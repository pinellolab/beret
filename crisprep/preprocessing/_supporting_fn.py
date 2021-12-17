from typing import List, Union, Literal
import subprocess as sb
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from .Edit import Allele, Edit

def _base_edit_to_from(start_base: chr = "A"):
    try:
        base_map = {"A":"G", "C":"T"}
    except KeyError:
        raise ValueError("Only A/C are supported for base to be edited.")
    return(base_map[start_base])


def _read_count_match(R1_filename: str, R2_filename: str) -> int:
    R1_count = _get_n_reads_fastq(R1_filename)
    R2_count = _get_n_reads_fastq(R2_filename)
    if R1_count != R2_count: 
        raise InputFileError("Paired end read numbers are different in R1({}) and R2({})".format(R1_count, R2_count)) 
    return(R1_count)


def _get_n_reads_fastq(fastq_filename: str):
    p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < %s | wc -l" % fastq_filename , shell=True,stdout=sb.PIPE)
    return(int(float(p.communicate()[0])/4.0))


def _get_fastq_handle(fastq_filename: str, mode = "r"):
    if fastq_filename.endswith('.gz'):
        fastq_handle=gzip.open(fastq_filename, mode)
    else:
        fastq_handle=open(fastq_filename, mode)
    return(fastq_handle)


def _read_is_good_quality(record: SeqIO.SeqRecord, 
        min_bp_quality = 0, 
        min_single_bp_quality = 0, 
        qend = -1):
    mean_quality_pass = np.array(record.letter_annotations["phred_quality"])[:qend_R1].mean() >= min_bp_quality
    min_quality_pass = np.array(record.letter_annotations["phred_quality"])[:qend_R1].min()>=min_single_bp_quality
    return(mean_quality_pass and min_quality_pass)


def _check_readname_match(R1:List[SeqIO.SeqRecord], R2:List[SeqIO.SeqRecord]):
   if len(R1) != len(R2):
    raise ValueError("The number of reads in R1 and R2 file does not match.")

    for i in range(len(R1)):
        R1_record = R1[i]
        R2_record = R2[i]

        if R1_record.name != R2_record.name : 
            raise InputFileError("R1 and R2 read discordance in read {} and {}".format(R1_record.name, R2_record.name))



def _get_guide_to_reporter_df(sgRNA_filename: str) -> pd.DataFrame:
    '''Returns a gRNA name to reporter sequence mapping.'''
    guide_to_reporter = {}

    with open(sgRNA_filename) as infile:
        sgRNA_df = pd.read_csv(infile)
        if not ('name' in sgRNA_df.columns and 'Reporter' in sgRNA_df.columns):
            raise InputFileError("Input gRNA file doesn't have the column 'gRNA' or 'gRNA_barcode'.")
        sgRNA_df.set_index('name', inplace = True)
        return(sgRNA_df)    


def revcomp(seq: Union[Seq, str]) -> str:
    if isinstance(seq, str): 
        seq = Seq(seq)
    return(str(seq.reverse_complement()))


def _fastq_iter_to_text(record):
    t, seq, q = record
    return("{}\n{}\n+\n{}\n".format(t, seq, q))


def _write_paired_end_reads(R1_record, R2_record, R1_out_handle, R2_out_handle):
        R1_out_handle.write(_fastq_iter_to_text(R1_record))
        R2_out_handle.write(_fastq_iter_to_text(R2_record))


def _get_edited_allele(
    ref_seq: str, 
    query_seq: str, 
    offset: int, 
    strand: Literal[1, -1] = 1,
    start_pos: int = 0, 
    end_pos: int = 100):

    allele = Allele()
    for i, (ref_nt, sample_nt) in enumerate(zip(ref_seq, query_seq)):
        if i < start_pos or i >= end_pos: continue
        if ref_nt == sample_nt: continue
        else: 
            edit = Edit(i - start_pos, ref_nt, sample_nt, offset, strand = strand)
            allele.add(edit)
    return(allele)



