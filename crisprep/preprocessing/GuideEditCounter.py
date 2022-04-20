from typing import Set, Literal, Iterable
import os
from os import path
import sys
import logging
from tqdm import tqdm
import pdb
import gzip
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqPhredIterator
from crisprep import ReporterScreen, Edit, Allele

# from crisprep.framework.Edit import Allele

from ._supporting_fn import (
    _base_edit_to_from,
    _read_count_match,
    _get_fastq_handle,
    _write_paired_end_reads,
    _get_edited_allele_crispresso,
    _check_readname_match,
    _read_is_good_quality,
    revcomp,
    _multiindex_dict_to_df,
    _write_alignment_matrix
)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n",
    datefmt="%a, %d %b %Y %H:%M:%S",
    stream=sys.stderr,
    filemode="w",
)
error = logging.critical
warn = logging.warning
debug = logging.debug
info = logging.info


class NTException(Exception):
    pass


class InputFileError(Exception):
    pass


class NoReadsAfterQualityFiltering(Exception):
    pass


class GuideEditCounter:
    def __init__(self, **kwargs):
        self.R1_filename = kwargs["R1"]
        self.R2_filename = kwargs["R2"]
        self.base_edited_from = kwargs["edited_base"]
        self.base_edited_to = _base_edit_to_from(self.base_edited_from)
        self.min_average_read_quality = kwargs["min_average_read_quality"]
        self.min_single_bp_quality = kwargs["min_single_bp_quality"]

        self.guides_info_df = pd.read_csv(kwargs["sgRNA_filename"])
        self.guides_has_strands = ("Strand" in self.guides_info_df.columns) or (
            "strand" in self.guides_info_df.columns
        )
        if self.guides_has_strands:
            info("Considering strand information of guides")
            assert "start_pos" in self.guides_info_df.columns
            # assert "gRNA_end" in self.guides_info_df.columns
        else:
            info("Ignoring guide strands, all guides are considered positive")

        self.qstart_R1 = kwargs["qstart_R1"]
        self.qend_R1 = kwargs["qend_R1"]
        self.qstart_R2 = kwargs["qstart_R2"]
        self.qend_R2 = kwargs["qend_R2"]
        self.name = kwargs["name"]

        self.sgRNA_filename = kwargs["sgRNA_filename"]
        self.count_only_bcmatched = False
        self._set_sgRNA_df()

        self.database_id = self._get_database_name()
        self.output_dir = os.path.join(
            os.path.abspath(kwargs["output_folder"]),
            "crisprep_count_%s" % self.database_id,
        )
        self._write_start_log()

        self.screen = ReporterScreen(
            X=np.zeros((len(self.guides_info_df), 1)),
            X_edit=np.zeros((len(self.guides_info_df), 1)),
            X_bcmatch=np.zeros((len(self.guides_info_df), 1)),
            guides=self.guides_info_df,
            condit=pd.DataFrame(index=[self.database_id]),
        )
        self.count_guide_edits = kwargs["count_guide_edits"]
        if self.count_guide_edits:
            self.screen.uns["guide_edit_counts"] = dict()
        self.count_reporter_edits = kwargs["count_reporter"]
        if self.count_reporter_edits:
            self.screen.uns["edit_counts"] = dict()
            self.gstart_reporter = kwargs["gstart_reporter"]

        self.count_edited_alleles = kwargs["count_allele"]
        if self.count_edited_alleles:
            self.screen.uns["allele_counts"] = dict()
            # pd.DataFrame(
            #     columns = ["guide_name", "allele_id", "count"]).set_index(
            #     ["guide_name", "allele_id"], drop = True)
            #self.screen.uns["alleles"] = pd.DataFrame(columns=["allele_id", "edit_id"])
            # Dict[guide_name -> List[edit_ID]]

        self.count_guide_reporter_alleles = kwargs["count_guide_reporter_alleles"]
        if self.count_guide_reporter_alleles:
            self.guide_to_guide_reporter_allele = dict()

        self.guide_start_seq = kwargs["guide_start_seq"]
        self.guide_bc = kwargs["guide_bc"]
        if self.guide_bc:
            self.guide_bc_len = kwargs["guide_bc_len"]

        self.offset = kwargs["offset"]
        if self.count_reporter_edits or self.count_edited_alleles or self.count_guide_reporter_alleles:
            self.reporter_length = kwargs["reporter_length"]
        self.guide_to_allele = dict()
        self.n_total_reads = _read_count_match(self.R1_filename, self.R2_filename)

        self.keep_intermediate = kwargs["keep_intermediate"]
        self.semimatch = 0
        self.bcmatch = 0
        self.nomatch = 0
        self.duplicate_match = 0
        self.duplicate_match_wo_barcode = 0

    def masked_equal(self, seq1, seq2):
        return(seq1.replace(self.base_edited_from, self.base_edited_to) == \
            seq2.replace(self.base_edited_from, self.base_edited_to))

    def _set_sgRNA_df(self):
        with open(self.sgRNA_filename) as infile:
            sgRNA_df = pd.read_csv(infile)
            if not ("name" in sgRNA_df.columns and "sequence" in sgRNA_df.columns):
                raise InputFileError(
                    "Input gRNA info file doesn't have the column 'name' or 'sequence'."
                )
            if self.count_only_bcmatched and "barcode" not in sgRNA_df.columns:
                raise InputFileError(
                    "Input gRNA info file doesn't have the column 'barcode'."
                )
            sgRNA_df = sgRNA_df.set_index("name")
        self.guides_info_df = sgRNA_df
        self.guides_info_df['masked_sequence'] = self.guides_info_df.sequence.map(lambda s: s.replace(self.base_edited_from, self.base_edited_to))
        self.guides_info_df['masked_barcode'] = self.guides_info_df.barcode.map(lambda s: s.replace(self.base_edited_from, self.base_edited_to))
        self.guide_lengths = sgRNA_df.sequence.map(lambda s: len(s)).unique()


    def check_filter_fastq(self):
        self.filtered_R1_filename = self._jp(
            os.path.basename(self.R1_filename).replace(".fastq", "").replace(".gz", "")
            + "_filtered.fastq.gz"
        )
        self.filtered_R2_filename = self._jp(
            os.path.basename(self.R2_filename).replace(".fastq", "").replace(".gz", "")
            + "_filtered.fastq.gz"
        )
        if path.exists(self.filtered_R1_filename) and path.exists(
            self.filtered_R2_filename
        ):
            warn("Using preexisting filtered file")
        else:
            self._check_names_filter_fastq()

    def get_counts(self):
        infile_R1 = _get_fastq_handle(self.R1_filename)
        infile_R2 = _get_fastq_handle(self.R2_filename)

        self.nomatch_R1_filename = self.R1_filename.replace(".fastq", "_nomatch.fastq")
        self.nomatch_R2_filename = self.R2_filename.replace(".fastq", "_nomatch.fastq")
        self.semimatch_R1_filename = self.R1_filename.replace(
            ".fastq", "_semimatch.fastq"
        )
        self.semimatch_R2_filename = self.R2_filename.replace(
            ".fastq", "_semimatch.fastq"
        )
        if self.count_edited_alleles or self.count_guide_reporter_alleles:
            _write_alignment_matrix(self.base_edited_from, self.base_edited_to, self.output_dir + "/.aln_mat.txt")

        if self.count_only_bcmatched:
            # count X
            self._get_guide_counts_bcmatch()
        else:  # count both bc matched & unmatched guide counts
            self._get_guide_counts_bcmatch_semimatch()

        if self.count_guide_edits:
            self.screen.uns["guide_edit_counts"] = _multiindex_dict_to_df(
                self.screen.uns["guide_edit_counts"], "edit", self.database_id
            )
            self.screen.uns["guide_edit_counts"].guide = self.screen.guides.index[
                self.screen.uns["guide_edit_counts"].guide.to_numpy(dtype=int)
            ]

        if self.count_reporter_edits:
            self.screen.uns["edit_counts"] = pd.DataFrame.from_dict(self.screen.uns["edit_counts"])

        if self.count_edited_alleles:
            guides = []
            alleles = []
            counts = []
            for guide, allele_to_count in self.guide_to_allele.items():
                if len(allele_to_count.keys()) == 0: continue
                guides.extend([guide]*len(allele_to_count.keys()))
                for allele, count in allele_to_count.items():
                    alleles.append(allele)
                    counts.append(count)
            if not (len(guides) == len(alleles) == len(counts)): 
                raise ValueError("Guides:{}, alleles:{}, counts:{}".format(len(guides), len(alleles), len(counts)))
            self.screen.uns["allele_counts"] = pd.DataFrame({"guide": guides, "allele": alleles, self.database_id: counts})
            
            if 'guide' in self.screen.uns["allele_counts"].columns:
                self.screen.uns["allele_counts"].guide = self.screen.guides.index[self.screen.uns["allele_counts"].guide]

        if self.count_guide_reporter_alleles:
            guides = []
            guide_alleles = []
            reporter_alleles = []
            counts = []
            for guide, allele_to_count in self.guide_to_guide_reporter_allele.items():
                if len(allele_to_count.keys()) == 0: continue
                guides.extend([guide]*len(allele_to_count.keys()))
                for key, count in allele_to_count.items():
                    reporter_allele, guide_allele = key
                    guide_alleles.append(guide_allele)
                    reporter_alleles.append(reporter_allele)
                    counts.append(count)
            if not (len(guides) == len(reporter_alleles) == len(guide_alleles) == len(counts)): 
                raise ValueError("Guides:{}, guide_alleles:{}, reporter_alleles: {}, counts:{}".format(
                    len(guides), len(guide_alleles), len(reporter_alleles), len(counts)))
            self.screen.uns["guide_reporter_allele_counts"] = pd.DataFrame(
                {"guide": guides, "guide_allele": guide_alleles, "reporter_allele": reporter_alleles, self.database_id: counts})
            
            if 'guide' in self.screen.uns["guide_reporter_allele_counts"].columns:
                self.screen.uns["guide_reporter_allele_counts"].guide = self.screen.guides.index[self.screen.uns["guide_reporter_allele_counts"].guide]            
        
        count_stat_path = self._jp("mapping_stats.txt")
        count_stat_file = open(count_stat_path, "w")
        count_stat_file.write("Read count with \n")
        count_stat_file.write("Unique guide match without barcode:\t{}\n".format(self.semimatch))
        count_stat_file.write("Unique guide match with barcode:\t{}\n".format(self.bcmatch))
        count_stat_file.write("No match:\t{}\n".format(self.nomatch))
        count_stat_file.write("Duplicate match with barcode:\t{}\n".format(self.duplicate_match))
        count_stat_file.write("No match with barcode & Duplicate match w/o barcode:\t{}\n".format(self.duplicate_match_wo_barcode))

    def _get_guide_counts_bcmatch(self):
        NotImplemented

    def _count_guide_edits(self, matched_guide_idx, R1_record: SeqIO.SeqRecord, single_base_qual_cutoff = 30):
        strand_str_to_int = {"neg": -1, "pos": 1}
        if self.guides_has_strands:
            try:
                guide_strand = strand_str_to_int[
                    self.screen.guides.Strand[matched_guide_idx]
                ]
            except KeyError:  # control guides
                guide_strand = 1
        else:
            guide_strand = 1
        ref_guide_seq = self.screen.guides.sequence[matched_guide_idx]
        read_guide_seq, read_guide_qual = self.get_guide_seq_qual(R1_record, len(ref_guide_seq))
        guide_edit_allele = read_guide_seq
        return(guide_edit_allele)
        

    def _count_reporter_edits(self, matched_guide_idx: int, R1_seq, R2_record: SeqIO.SeqRecord, 
    single_base_qual_cutoff=30, guide_allele = None,):
        '''
        Count edits in single read to save as allele.
        '''
        strand_str_to_int = {"neg": -1, "pos": 1}
        ref_reporter_seq = self.screen.guides.Reporter[matched_guide_idx]
        read_reporter_seq, read_reporter_qual = self.get_reporter_seq_qual(R2_record)

        if self.guides_has_strands:
            try:
                guide_strand = strand_str_to_int[
                    self.screen.guides.Strand[matched_guide_idx]
                ]
                if guide_strand == -1:
                    offset = (
                        self.screen.guides.start_pos[matched_guide_idx]
                        + self.screen.guides.guide_len[matched_guide_idx]
                        - 1
                    )
                if guide_strand == 1:
                    offset = self.screen.guides.start_pos[matched_guide_idx]
            except KeyError:  # control guides
                guide_strand = 1
                offset = 0
        else:
            guide_strand = 1
            offset = -(self.screen.guides["Target base position in reporter"][matched_guide_idx] -1)
            # TODO: clean this up
        
        allele = read_reporter_seq
        print(allele)
        if True:
            if self.count_edited_alleles:
                if matched_guide_idx in self.guide_to_allele.keys():
                    if allele in self.guide_to_allele[matched_guide_idx].keys():
                        self.guide_to_allele[matched_guide_idx][allele] += 1
                    else:
                        self.guide_to_allele[matched_guide_idx][allele] = 1
                else:
                    self.guide_to_allele[matched_guide_idx] = {allele: 1}
                #self._write_allele(matched_guide_idx, allele)
            if self.count_guide_reporter_alleles and (not guide_allele is None):
                if matched_guide_idx in self.guide_to_guide_reporter_allele.keys():
                    if (allele, guide_allele) in self.guide_to_guide_reporter_allele[matched_guide_idx].keys():
                        self.guide_to_guide_reporter_allele[matched_guide_idx][(allele, guide_allele)] += 1
                    else:
                        self.guide_to_guide_reporter_allele[matched_guide_idx][(allele, guide_allele)] = 1
                else:
                    self.guide_to_guide_reporter_allele[matched_guide_idx] = {(allele, guide_allele): 1}
        #self._write_edits(matched_guide_idx, allele)

    def _get_guide_counts_bcmatch_semimatch(
        self, bcmatch_layer="X_bcmatch", semimatch_layer="X"
    ):

        self.screen.layers[semimatch_layer] = np.zeros_like((self.screen.X))
        R1_iter, R2_iter = self._get_fastq_iterators()

        outfile_R1_nomatch, outfile_R2_nomatch = self._get_fastq_handle("nomatch")
        outfile_R1_semimatch, outfile_R2_semimatch = self._get_fastq_handle("semimatch")
        outfile_R1_dup_wo_bc, outfile_R2_dup_wo_bc = self._get_fastq_handle(
            "duplicate_wo_barcode"
        )
        outfile_R1_dup, outfile_R2_dup = self._get_fastq_handle("duplicate")

        for i, (r1, r2) in tqdm(
            enumerate(zip(R1_iter, R2_iter)), total=self.n_reads_after_filtering
        ):
            R1_seq = str(r1.seq)
            R2_seq = str(r2.seq)

            bc_match, semimatch = self._match_read_to_sgRNA_bcmatch_semimatch(
                R1_seq, R2_seq
            )

            if len(bc_match) == 0:
                if len(semimatch) == 0:  # no guide match
                    if self.keep_intermediate:
                        _write_paired_end_reads(
                            r1, r2, outfile_R1_nomatch, outfile_R2_nomatch
                        )
                    self.nomatch += 1
                elif len(semimatch) >= 2:  # Duplicate match if w/o barcode
                    if self.keep_intermediate:
                        _write_paired_end_reads(
                            r1, r2, outfile_R1_dup_wo_bc, outfile_R2_dup_wo_bc
                        )
                    self.duplicate_match_wo_barcode += 1
                else:  # guide match with no barcode match
                    matched_guide_idx = semimatch[0]
                    try:
                        self.screen.layers[semimatch_layer][matched_guide_idx, 0] += 1
                    except:
                        print(semimatch)

                    if self.count_guide_edits:
                        self._count_guide_edits(matched_guide_idx, r1)
                    self.semimatch += 1

            elif len(bc_match) >= 2:  # duplicate mapping
                if self.keep_intermediate:
                    _write_paired_end_reads(r1, r2, outfile_R1_dup, outfile_R2_dup)
                self.duplicate_match += 1

            else:  # unique barcode match
                matched_guide_idx = bc_match[0]
                self.screen.layers[bcmatch_layer][matched_guide_idx, 0] += 1
                self.bcmatch += 1
                if self.count_guide_edits or self.count_guide_reporter_alleles:
                    guide_allele = self._count_guide_edits(matched_guide_idx, r1)
                if self.count_reporter_edits or self.count_edited_alleles or self.count_guide_reporter_alleles:
                    # TBD: what if reporter seq doesn't match barcode & guide?
                    if self.count_guide_reporter_alleles:
                        self._count_reporter_edits(matched_guide_idx, R1_seq, r2, guide_allele = guide_allele)
                    else: 
                        self._count_reporter_edits(matched_guide_idx, R1_seq, r2)

        self.screen.X = (
            self.screen.layers[semimatch_layer] + self.screen.layers[bcmatch_layer]
        )

    def _write_allele(self, guide_idx: int, allele: Allele, uns_key="allele_counts"):
        # if len(allele.edits) == 0:
        #     return

        if (guide_idx, str(allele)) in self.screen.uns[uns_key].keys():
            self.screen.uns[uns_key][(guide_idx, str(allele))] += 1
        else:
            self.screen.uns[uns_key][(guide_idx, str(allele))] = 1

    def _write_edits(self, guide_idx: int, allele: Allele, uns_key="edit_counts"):
        for edit in allele.edits:
            if (guide_idx, str(edit)) in self.screen.uns[uns_key].keys():
                self.screen.uns[uns_key][(guide_idx, str(edit))] += 1
            else:
                self.screen.uns[uns_key][(guide_idx, str(edit))] = 1

    def _write_guide_reporter_allele(self, guide_idx: int, allele: Allele, guide_allele: Allele, uns_key="guide_reporter_allele_counts"):
        if len(allele.edits) == 0 and len(guide_allele.edits) == 0:
            return
        if (guide_idx, str(allele), str(guide_allele)) in self.screen.uns[uns_key].keys():
            self.screen.uns[uns_key][(guide_idx, str(allele), str(guide_allele))] += 1
        else:
            self.screen.uns[uns_key][(guide_idx, str(allele), str(guide_allele))] = 1

    def get_guide_seq(self, R1_seq, R2_seq, guide_length):
        """This can be edited by user based on the read construct."""
        guide_start_idx = R1_seq.find(self.guide_start_seq)
        if guide_start_idx == -1:
            return None
        if guide_start_idx + guide_length >= len(R1_seq): 
            return None
        guide_start_idx = guide_start_idx + len(self.guide_start_seq)
        gRNA_seq = R1_seq[guide_start_idx : (guide_start_idx + guide_length)]
        if len(gRNA_seq) != guide_length:
            return None
        return gRNA_seq

    def get_guide_seq_qual(self, R1_record: SeqIO.SeqRecord, guide_length):
        R1_seq = R1_record.seq
        guide_start_idx = R1_seq.find(self.guide_start_seq)
        if guide_start_idx == -1:
            return None, None
        if guide_start_idx + guide_length >= len(R1_seq): 
            return None, None
        guide_start_idx = guide_start_idx + len(self.guide_start_seq)
        seq = R1_record[guide_start_idx : (guide_start_idx + guide_length)]
        return(str(seq.seq), seq.letter_annotations["phred_quality"])

    def get_reporter_seq(self, R1_seq, R2_seq):
        """This can be edited by user based on the read construct."""
        reporter_seq = revcomp(
            R2_seq[self.guide_bc_len : (self.guide_bc_len + self.reporter_length)]
        )
        return(reporter_seq)

    def get_reporter_seq_qual(self, R2_record: SeqIO.SeqRecord):
        seq = R2_record[self.guide_bc_len : (self.guide_bc_len + self.reporter_length)].reverse_complement()
        return(str(seq.seq), seq.letter_annotations["phred_quality"])

    def get_barcode(self, R1_seq, R2_seq):
        """This can be edited by user based on the read construct."""
        return revcomp(R2_seq[: self.guide_bc_len])

    def _match_read_to_sgRNA_bcmatch_semimatch(self, R1_seq: str, R2_seq: str):
        # This should be adjusted for each experimental recipes.'
        guide_barcode = self.get_barcode(R1_seq, R2_seq)
        bc_match_idx = np.array([])
        semimatch_idx = np.array([])
        for guide_length in self.guide_lengths:
            seq = self.get_guide_seq(R1_seq, R2_seq, guide_length)
            if seq is None:
                continue

            _seq_match = np.where(seq.replace(self.base_edited_from, self.base_edited_to) == self.screen.guides.masked_sequence)[0]
            _bc_match = np.where(guide_barcode.replace(self.base_edited_from, self.base_edited_to) == self.screen.guides.masked_barcode)[0]
            
            bc_match_idx = np.append(bc_match_idx, np.intersect1d(_seq_match, _bc_match))
            semimatch_idx = np.append(semimatch_idx, _seq_match)

        return (bc_match_idx.astype(int), semimatch_idx.astype(int))

    def _get_guide_position_seq_of_read(self, seq):
        guide_start_idx = self._get_guide_start_idx(seq)
        if guide_start_idx == -1:
            return None

        return [
            seq[guide_start_idx : (guide_start_idx + guide_length)]
            for guide_length in self.guide_lengths
        ]

    def _get_guide_start_idx(self, seq):
        start_seq_idx = seq.find(self.guide_start_seq)
        if start_seq_idx == -1:
            return -1
        return start_seq_idx + len(self.guide_start_seq)

    def get_gRNA_barcode(self, R1_seq, R2_seq):
        # This can be adjusted for different construct design.
        return revcomp(R2_seq[: self.guide_bc_len])

    def _get_fastq_handle(
        self,
        out_type: Literal[
            "semimatch", "nomatch", "duplicate_wo_barcode", "duplicate"
        ] = None,
    ):
        R1_filename = self.R1_filename.replace(".fastq", "_{}.fastq".format(out_type))
        R2_filename = self.R2_filename.replace(".fastq", "_{}.fastq".format(out_type))
        R1_handle = _get_fastq_handle(R1_filename, "w")
        R2_handle = _get_fastq_handle(R2_filename, "w")

        return (R1_handle, R2_handle)

    def _get_fastq_iterators(self):
        R1_handle = _get_fastq_handle(self.R1_filename)
        R2_handle = _get_fastq_handle(self.R2_filename)

        R1_iterator = FastqPhredIterator(R1_handle)
        R2_iterator = FastqPhredIterator(R2_handle)

        return (R1_iterator, R2_iterator)

    def _get_seq_records(self):
        R1_handle = _get_fastq_handle(self.R1_filename)
        R2_handle = _get_fastq_handle(self.R2_filename)
        R1 = list(SeqIO.parse(R1_handle, "fastq"))
        R2 = list(SeqIO.parse(R2_handle, "fastq"))
        R1_handle.close()
        R2_handle.close()
        return (R1, R2)

    def _check_names_filter_fastq(self, filter_by_qual=False):
        if self.min_average_read_quality > 0 or self.min_single_bp_quality > 0:
            info(
                "Filtering reads with average bp quality < {} and single bp quality < {} ...".format(
                    self.min_average_read_quality, self.min_single_bp_quality
                )
            )

        if self.qend_R1 > 0 or self.qend_R2 > 0:
            info(
                "In the filtering, bases up to position {} of R1 and {} of R2 are considered.".format(
                    self.qend_R1, self.qend_R2
                )
            )

        R1, R2 = self._get_seq_records()

        _check_readname_match(R1, R2)
        if filter_by_qual:
            self.n_reads_after_filtering = self._filter_read_quality(R1, R2)
            if self.n_reads_after_filtering == 0:
                raise NoReadsAfterQualityFiltering(
                    "No reads in input or no reads survived the average or single bp quality filtering."
                )
            else:
                info(
                    "Number of reads in input:%d\tNumber of reads after filtering:%d\n"
                    % (self.n_total_reads, self.n_reads_after_filtering)
                )
        else:
            self.n_reads_after_filtering = self.n_total_reads

    def _filter_read_quality(self, R1=None, R2=None) -> int:
        R1_filtered = gzip.open(self.filtered_R1_filename, "w+")
        R2_filtered = gzip.open(self.filtered_R2_filename, "w+")

        if R1 is None or R2 is None:
            R1, R2 = self._get_seq_records()

        n_reads_after_filtering = 0
        for i in range(len(R1)):
            R1_record = R1[i]
            R2_record = R2[i]

            if self.filter_by_qual:
                R1_quality_pass = _read_is_good_quality(
                    R1_record,
                    self.min_average_read_quality,
                    self.min_single_bp_quality,
                    self.qend_R1,
                )
                R2_quality_pass = _read_is_good_quality(
                    R2_record,
                    self.min_average_read_quality,
                    self.min_single_bp_quality,
                    self.qend_R2,
                )

                if R1_quality_pass and R2_quality_pass:
                    n_reads_after_filtering += 1
                    R1_filtered.write(R1.format("fastq"))
                    R2_filtered.write(R2.format("fastq"))
        return n_reads_after_filtering

    def _write_start_log(self):
        try:
            os.makedirs(self.output_dir)
            print("Creating Folder %s" % self.output_dir)
        except:
            print("Folder %s already exists." % self.output_dir)
        self.log_filename = self._jp("CRISPRepCount_RUNNING_LOG.txt")
        logging.getLogger().addHandler(logging.FileHandler(self.log_filename))
        with open(self.log_filename, "w+") as outfile:
            outfile.write(
                "[Command used]:\nCRISPRessoCount %s\n\n[Execution log]:\n"
                % " ".join(sys.argv)
            )

    def _jp(self, filename):
        return os.path.join(self.output_dir, filename)

    def _get_database_name(self):
        get_name_from_fasta = (
            lambda x: os.path.basename(x)
            .replace(".fastq", "")
            .replace(".gz", "")
            .replace("_R1", "")
        )

        if not self.name:
            database_id = "%s" % get_name_from_fasta(self.R1_filename)
        else:
            database_id = self.name
        return database_id
