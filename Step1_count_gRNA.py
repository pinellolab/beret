# -*- coding: utf-8 -*-

import os
from os import path
import gzip
import argparse
import sys
import gzip
import subprocess as sb
from collections import defaultdict
import unicodedata
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from typing import Dict, List, Union, Tuple
from tqdm import tqdm 

import logging
logging.basicConfig(level=logging.INFO,
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info


_ROOT = os.path.abspath(os.path.dirname(__file__))


####Support functions###

def check_file(filename):
    try:
        with open(filename): pass
    except IOError:
        raise Exception('I cannot open the file: '+filename)
 

def check_library(library_name):
        try:
                return __import__(library_name)
        except:
                error('You need to install %s module to use CRISPRessoCount!' % library_name)
                sys.exit(1)


def slugify(value): #adapted from the Django project
    
    value = unicodedata.normalize('NFKD', unicode(value)).encode('ascii', 'ignore')
    value = unicode(re.sub('[^\w\s-]', '_', value).strip())
    value = unicode(re.sub('[-\s]+', '-', value))
    
    return str(value)


def get_fastq_handle(fastq_filename):
    if fastq_filename.endswith('.gz'):
            fastq_handle=gzip.open(fastq_filename)
    else:
            fastq_handle=open(fastq_filename)
    return(fastq_handle)


def read_is_good_quality(record: SeqIO.SeqRecord, 
        min_bp_quality = 0, 
        min_single_bp_quality = 0, 
        qend = -1):
    mean_quality_pass = np.array(record.letter_annotations["phred_quality"])[:qend_R1].mean() >= min_bp_quality
    min_quality_pass = np.array(record.letter_annotations["phred_quality"])[:qend_R1].min()>=min_single_bp_quality
    return(mean_quality_pass and min_quality_pass)


def check_read_names_and_filter_quality(R1_filename, 
        R2_filename, 
        output_R1_filename=None, 
        output_R2_filename = None,
        filter_by_qual = False,
        min_bp_quality=20,
        min_single_bp_quality=0, 
        qend_R1 = -1, 
        qend_R2 = -1):
        
        if min_bp_quality > 0 or min_single_bp_quality > 0:
            info('Filtering reads with average bp quality < {} \
                and single bp quality < {} ... in up to position {} of R1 and {} of R2'.format(args.min_average_read_quality,
                    args.min_single_bp_quality))

        if qend_R1 > 0 or qend_R2 > 0:
            info("In the filtering, bases up to position {} of R1 and {} of R2 are considered.".format(qend_R1, qend_R2))

        R1_handle = get_fastq_handle(R1_filename)
        R2_handle = get_fastq_handle(R2_filename)

        if not output_R1_filename:
                output_R1_filename = R1_filename.replace('.fastq', '').replace('.gz', '') + '_filtered.fastq.gz'
                output_R2_filename = R2_filename.replace('.fastq', '').replace('.gz', '') + '_filtered.fastq.gz'

        try:
            if filter_by_qual:
                R1_filtered = gzip.open(output_R1_filename, 'w+')
                R2_filtered = gzip.open(output_R2_filename, 'w+')
            R1 = list(SeqIO.parse(R1_handle, "fastq"))
            R2 = list(SeqIO.parse(R2_handle, "fastq"))
            R1_handle.close()
            R2_handle.close()

            if len(R1) != len(R2):
                raise ValueError("The number of reads in R1 and R2 file does not match.")

            for i in range(len(R1)):
                R1_record = R1[i]
                R2_record = R2[i]

                if R1_record.name != R2_record.name : 
                    raise InputFileError("R1 and R2 read discordance in read {} and {}".format(R1_record.name, R2_record.name))
                if filter_by_qual:
                    R1_quality_pass = read_is_good_quality(R1_record, min_bp_quality, min_single_bp_quality, qend_R1)
                    R2_quality_pass = read_is_good_quality(R2_record, min_bp_quality, min_single_bp_quality, qend_R2)

                    if R1_quality_pass and R2_quality_pass:
                        R1_filtered.write(record.format('fastq'))
                        R2_filtered.write(record.format('fastq'))
        except:
                raise Exception('Error handling the fastq_filtered_outfile')

        if filter_by_qual: 
            return((output_R1_filename, output_R2_filename))
        else:
            return((R1_filename, R2_filename))


def find_wrong_nt(sequence):
    return(list(set(sequence.upper()).difference(set(['A','T','C','G','N']))))


def get_n_reads_fastq(fastq_filename):
     p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < %s | wc -l" % fastq_filename , shell=True,stdout=sb.PIPE)
     return(int(float(p.communicate()[0])/4.0))


def mask_sequence_positions(seq: str, pos: np.array) -> str:
    return("".join([seq[i] if i in pos else "N" for i in range(len(seq))]))


def revcomp(seq: Union[Seq, str]) -> str:
    if isinstance(seq, str): 
        seq = Seq(seq)
    return(str(seq.reverse_complement()))


###EXCEPTIONS############################
class NTException(Exception):
    pass

class InputFileError(Exception):
    pass

class NoReadsAfterQualityFiltering(Exception):
    pass
#########################################

def get_input_parser():
    """Get the input data"""
    print('  \n~~~CRISPRessoCount~~~')
    print('-Utility to perform sgRNA and reporter count from CRISPR base editors-')
    print(r'''
          )                                             )
         (           ________________________          (
        __)__       | __   __            ___ |        __)__
     C\|     \      |/  ` /  \ |  | |\ |  |  |     C\|     \
       \     /      |\__, \__/ \__/ | \|  |  |       \     /
        \___/       |________________________|        \___/
    ''')
    
    
    print('\n[Luca Pinello 2017, Jayoung Ryu 2021, send bugs, suggestions or *green coffee* to jayoung_ryu AT g DOT harvard DOT edu]\n\n')
    
    
    parser = argparse.ArgumentParser(description='CRISPRessoCount parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--R1', type=str,  help='fastq file for read 1', required=True,default='Fastq filename' )
    parser.add_argument('--R2', type=str,  help='fastq file for read 2, sorted as the same name order as in --R1 file.', required=True, default='Fastq filename' )
    parser.add_argument('-b', '--edited_base', type = str, required = True, help = 'For base editors, the base that should be ignored when matching the gRNA sequence')
    parser.add_argument('-f','--sgRNA_file', type=str, required = True, help='''sgRNA description file. The format requires three columns: gRNA, Reporter, gRNA_barcode.''')

    #optional
    parser.add_argument('--guide_start', type = str, help = "Guide starts after this sequence in R1", default = "GGAAAGGACGAAACACCG")
    parser.add_argument('-r', '--count_reporter', help = "Count reporter edits.", action = 'store_true')
    parser.add_argument('-q','--min_average_read_quality', type=int, help='Minimum average quality score (phred33) to keep a read', default=0)
    parser.add_argument('-s','--min_single_bp_quality', type=int, help='Minimum single bp score (phred33) to keep a read', default=0)
    parser.add_argument('-n','--name',  help='Output name', default='')
    parser.add_argument('-o','--output_folder',  help='', default='')
    parser.add_argument('-l', '--reporter_length', type = int, help = "length of the reporter", default = 32)
    parser.add_argument('--keep_intermediate',help='Keep all the  intermediate files',action='store_true')
    parser.add_argument('--qstart_R1', help='Start position of the read when filtering for quality score of the read 1', type = int, default = 0)
    parser.add_argument('--qend_R1', help = 'End position of the read when filtering for quality score of the read 1', type = int, default = 47)
    parser.add_argument('--qstart_R2', help = 'Same as qstart_R1, for read 2 fastq file', default = 0)
    parser.add_argument('--qend_R2', help = 'Same as qstart_R2, for read 2 fastq file', default = 36)
    parser.add_argument('--guide_bc_len', help = 'Guide barcode sequence length at the beginning of the R2', type = str, default = 4)
    parser.add_argument('--offset', help = 'Guide file has offest column that will be added to the relative position of reporters.', action = 'store_true')
    parser.add_argument('--align_fasta', help = 'gRNA is aligned to this sequence to infer the offset. Can be used when the exact offset is not provided.', type = str, default = '')

    return(parser)


def check_arguments(args):
    """Check the argument validity of the ArgumentParser"""
    check_file(args.R1)
    check_file(args.R2)

    if args.sgRNA_file:
        check_file(args.sgRNA_file)

    # Edited base should be one of A/C/T/G
    if args.edited_base.upper() not in ['A', 'C', "T", "G"]:
        raise ValueError("The edited base should be one of A/C/T/G, {} provided.".format(args.edited_base))
    else:
        edited_base = args.edited_base.upper()
    info('Using specified edited base: %s' % edited_base)

    read_length = get_first_read_length(args.R1)

    # Check if positions of guide and quality control is valid
    NotImplemented

    if not (args.qstart_R1 < read_length and args.qstart_R2 < read_length):
        raise ValueError("The start position of base quality filter is not nonnegative ({} for R1, {} for R2 provided)".format(args.qstart_R1, args.qstart_R2))
    
    if not (args.qend_R1 < read_length and args.qend_R2 < read_length):
        raise ValueError("The start position of base quality filter is not nonnegative ({} for R1, {} for R2 provided)".format(args.qstart_R1, args.qstart_R2))

    if args.qend_R2 != args.guide_bc_len + args.reporter_length :
        warn("Quality of R2 checked up until {}bp, while the length of guide barcode and reporter combined is {}bp.".format(
            args.qend_R2, args.guide_bc_len + args.reporter_length))
    info("Using guide barcode length {}, guide start {}, and guide length {}".format(args.guide_bc_len, args.guide_start, ','.join(list(map(lambda x: str(x), GUIDE_LENGTHS))))) 
    #normalize name and remove not allowed characters
    if args.name:   
        clean_name=slugify(args.name)
        if args.name!= clean_name:
               warn('The specified name %s contained characters not allowed and was changed to: %s' % (args.name,clean_name))
               args.name=clean_name

    if args.offset:
        df = pd.read_csv(args.sgRNA_file)
        if not 'offset' in df.columns:
            raise InputFileError("Offset option is set but the input file doesn't contain the offset column.")
        if len(args.align_fasta) > 0:
            error("Can't have --offset and --align_fasta option together.")

    info("Done checking input arguments.")

    return(args)

    #count reads
#def check_input_reads(args):
#    N_READS_INPUT=get_n_reads_fastq(args.R1)
#    N_READS_INPUT_R2 = get_n_reads_fastq(args.R2)
#    if N_READS_INPUT != N_READS_INPUT_R2 : 
#        raise InputFileError("Number of reads not matched in two input fastq file")
#    return(N_READS_INPUT)



def get_database_id(args):    
    get_name_from_fasta = lambda  x: os.path.basename(x).replace('.fastq','').replace('.gz','').replace("_R1",'')
    
    if not args.name:
        database_id='%s' % get_name_from_fasta(args.R1)    
    else:
        database_id=args.name
    return(database_id)


def read_count_match(R1_filename, R2_filename) -> int:
    R1_count = get_n_reads_fastq(R1_filename)
    R2_count = get_n_reads_fastq(R2_filename)
    if R1_count != R2_count: 
        raise InputFileError("Paired end read numbers are different in R1({}) and R2({})".format(R1_count, R2_count)) 
    return(R1_count)


def get_first_read_length(fastq_filename):
    for record in SeqIO.parse(fastq_filename, "fastq"):
        return(len(record))
    raise InputFileError("Provided R1 file doesn't have any read to parse")


def get_sgRNA_dict(sgRNA_filename: str, edited_base: str) -> Dict[str, Dict[np.array, Dict[str, str]]]:
    """Reads the sgRNA file (csv) with three columns: name (of gRNA), gRNA(sequence), gRNA_barcode.
    Will create a dictionary of (barcode -> (masked position -> (masked sequence -> gRNA name))) 
    for fast search of gRNA match with masked positon.
    """

    info('Using guides information in %s' % sgRNA_filename)
    masked_guides_dict = dict()
    guides = dict()

    with open(sgRNA_filename) as infile:
        sgRNA_df = pd.read_csv(infile)
        if not ('gRNA' in sgRNA_df.columns and 'gRNA_barcode' in sgRNA_df.columns):
            raise InputFileError("Input gRNA file doesn't have the column 'gRNA' or 'gRNA_barcode'.")
        
        duplicate_gRNA = 0
        for i, guide in enumerate(sgRNA_df["gRNA"]):
            if not len(guide) in GUIDE_LENGTHS:
                warn("Guide length does not match the prespecified length ({}bp)".format(",".join(GUIDE_LENGTHS)))
            guide = guide.replace(edited_base, "N")
            
            # Writing gRNA stats with duplication, masked gRNA
            if guide in guides.keys(): 
                guides[guide] += 1
            else: guides[guide] = 1

            gRNA_bc = sgRNA_df["gRNA_barcode"][i]
            gRNA_seq = sgRNA_df["gRNA"][i]
            mask_pos = tuple(np.nonzero(gRNA_seq == args.edited_base)[0])  # Indices of the masked position (np.array)
            masked_seq = gRNA_seq.replace(edited_base, "N")
            gRNA_name = sgRNA_df["name"][i]

            if gRNA_bc in masked_guides_dict.keys():
                if mask_pos in masked_guides_dict[gRNA_bc].keys():
                    if masked_seq in masked_guides_dict[gRNA_bc][mask_pos].keys():
                        duplicate_gRNA += 1
                    else:
                        masked_guides_dict[gRNA_bc][mask_pos][masked_seq] = gRNA_name
                else:
                    masked_guides_dict[gRNA_bc][mask_pos] = {masked_seq: gRNA_name}                    
            else: 
                masked_guides_dict[gRNA_bc] = {mask_pos : {masked_seq: gRNA_name}}
    return(masked_guides_dict)


def get_sgRNA_dict_simple(sgRNA_filename):
    seq_to_name = dict()
    # TBD: change this into dataframe
    with open(sgRNA_filename) as infile:
        sgRNA_df = pd.read_csv(infile)
        if not ('name' in sgRNA_df.columns and 'gRNA_barcode' in sgRNA_df.columns and 'gRNA' in sgRNA_df.columns):
            raise InputFileError("Input gRNA file doesn't have the column 'gRNA' or 'gRNA_barcode or 'name'.")
        for i in range(len(sgRNA_df)):
            seq_to_name[sgRNA_df["gRNA"][i]] = (sgRNA_df["name"][i], sgRNA_df["gRNA_barcode"][i])
    return(seq_to_name)


def get_guide_info(sgRNA_filename: str) -> pd.DataFrame:
    '''Returns a gRNA name to reporter sequence mapping.'''
    guide_to_reporter = {}

    with open(sgRNA_filename) as infile:
        sgRNA_df = pd.read_csv(infile)
        if not ('name' in sgRNA_df.columns and 'Reporter' in sgRNA_df.columns):
            raise InputFileError("Input gRNA file doesn't have the column 'gRNA' or 'gRNA_barcode'.")
        sgRNA_df.set_index('name', inplace = True)
        return(sgRNA_df)    
    

def match_masked_sgRNA(masked_guides_dict: dict, 
    gRNA_seq: str, 
    gRNA_barcode: str) -> List[str]:
    ''' 
    Match a gRNA_seq to the masked_guides_dict

    masked_guides_dict : dict
        (barcode -> (masked position -> (masked sequence -> gRNA name)))
    '''
    n_matches = 0
    matches = []
    if gRNA_barcode in masked_guides_dict.keys():
        for mask, maseq2name in masked_guides_dict[gRNA_barcode].items():   
            masked_seq = "".join(["N" if i in mask else gRNA_seq[i] for i in range(len(gRNA_seq))])
            # maseq2name : (masked sequence -> gRNA name)
            if masked_seq in maseq2name.keys():
                n_matches += 1
                matches.append(maseq2name[masked_seq])
            else:
                pass # no match with current mask. Try other mask
    else:
        pass # no barcode match.
    return(matches)


def gRNA_eq(guide, observed, edit_start, edit_end):
    if len(observed) != len(guide): return(False)
    for observed_nt, guide_nt in zip(observed, guide):
        if (observed_nt == guide_nt) or ((guide_nt, observed_nt) == (edit_start, edit_end)):
            continue
        else:
            return(False)
    return(True)


def match_sgRNA(guides_dict: Dict[str, str], query_seq: str, guide_bc : str, 
    edit_start: str = "A", edit_end: str = "G") -> List[str]:
    '''
    Search across all the gRNAs to find the match
    guides_dict: seq to name
    '''
    # TBD: condition check for edit_start, end
    matches = []
    for ref_seq, (name, ref_barcode) in guides_dict.items():
        if ref_barcode == guide_bc and gRNA_eq(ref_seq, query_seq, edit_start, edit_end) : 
#        if gRNA_eq(ref_seq, query_seq, edit_start, edit_end): 
            matches.append(name)
    return(matches)


def count_masked_guides(R1_filename, R2_filename, 
    guides_dict: Dict[str, str], 
    guide_start: str, 
    guide_bc_len: int, 
    write_nomatch: bool = False, 
    count_reporter_edits: bool = False, 
    guide_info_df: pd.DataFrame = None) -> Union[Dict[str, int], Dict[str, Dict[tuple, int]]]:
    '''
    Given a read pair, find matches among gRNA and optionally find the reporter edits.
    Returns a dictionary of (guide name -> count) if count_reporter_edits = False.
    Returns a dictionary of (guide name -> ((pos, ref_base, edited_base) -> count)) otherwise.
    '''

    infile_R1 = get_fastq_handle(R1_filename)
    infile_R2 = get_fastq_handle(R2_filename)

    if write_nomatch:
        nomatch_R1_filename = R1_filename.replace(".fastq", "_nomatch.fastq")
        nomatch_R2_filename = R2_filename.replace(".fastq", "_nomatch.fastq")
        outfile_R1_nomatch = open(nomatch_R1_filename, "w")
        outfile_R2_nomatch = open(nomatch_R2_filename, "w")

    if count_reporter_edits:
        assert not guide_info_df is None, "No guide to reporter dictionary passed onto count_masked_guides."

    info('Counting sgRNAs...')
    
    iterator_R1 = FastqGeneralIterator(infile_R1)
    iterator_R2 = FastqGeneralIterator(infile_R2)

    gname_to_count = dict()

    for i, (r1, r2) in tqdm(enumerate(zip(iterator_R1, iterator_R2)), total = N_READS_AFTER_PREPROCESSING):
#    for i, (r1, r2) in enumerate(zip(iterator_R1, iterator_R2)):
        t1, R1_seq, q1 = r1
        t2, R2_seq, q2 = r2
        gRNA_names = []
        for guide_length in GUIDE_LENGTHS:
            gRNA_barcode = revcomp(R2_seq[:guide_bc_len])
            guide_start_idx = R1_seq.find(guide_start)
            if guide_start_idx == -1 : continue
            guide_start_idx = guide_start_idx + len(guide_start)
            gRNA_seq = R1_seq[guide_start_idx:guide_start_idx + guide_length]
            gRNA_names.extend(match_sgRNA(guides_dict, gRNA_seq, gRNA_barcode))

        if len(gRNA_names) == 0: 
            if write_nomatch:
                outfile_R1_nomatch.write("{}\n{}\n+\n{}\n".format(t1, R1_seq, q1))
                outfile_R2_nomatch.write("{}\n{}\n+\n{}\n".format(t2, R2_seq, q2))
        elif len(gRNA_names) >= 2: 
            warn("{}th read with {} nonunique mapping to {}. Discarding the read.".format(i, len(gRNA_names), ",".join(gRNA_names)))
            print("{}\n{}\n+\n{}\n".format(t1, R1_seq, q1))
            print("{}\n{}\n+\n{}\n".format(t2, R2_seq, q2))
        else:
            # unique match
            gRNA_name = gRNA_names[0]

            if count_reporter_edits:
                if gRNA_name in gname_to_count.keys():
                    reporter_edit_counts = gname_to_count[gRNA_name]
                else: 
                    reporter_edit_counts = dict()
                
                reporter_ref = guide_info_df.loc[gRNA_name].Reporter
                reporter_seq = revcomp(R2_seq[guide_bc_len:guide_bc_len + REPORTER_LENGTH])
                reporter_seq_qual = q2[guide_bc_len:guide_bc_len + REPORTER_LENGTH][::-1]
                
                if args.offset :
                    offset = guide_info_df.loc[gRNA_name].offset
                
                elif len(args.align_fasta) > 0 :
                    try:
                        gene_seq = next(SeqIO.parse(args.align_fasta, "fasta")).seq.upper()
                        #warn("Treating the first entry of the file as the align reference")
                    except:
                        raise InputFileError("Please check the gene sequence fasta file.")
                    
                    guide_info_df["pos_gRNA_seq"] = guide_info_df.gRNA

                    if "Strand" in guide_info_df.columns:
                        guide_info_df.pos_gRNA_seq.loc[(guide_info_df.Strand == "neg") | (guide_info_df.Strand == '-')] = guide_info_df.pos_gRNA_seq.loc[guide_info_df.Strand == "neg"].apply(revcomp)
                    offset = gene_seq.find(guide_info_df.pos_gRNA_seq.loc[gRNA_name])
                
                # Check if gRNA can be mapped to multiple locations of the gene sequence
                    if gene_seq.count(guide_info_df.loc[gRNA_name].pos_gRNA_seq) > 1:
                        warn("gRNA {} can be mapped to multiple region. Using the first mapped position {}. The provided position is {}.".format(gRNA_name, offset, guide_info_df.loc[gRNA_name].pos))
                        #exit(1)
                    
                for i, (ref_nt, sample_nt) in enumerate(zip(reporter_ref, reporter_seq)):
                    if i < EDIT_START_POS or i >= EDIT_END_POS: continue
                    if ref_nt == sample_nt: continue
                    if ord(reporter_seq_qual[i]) - 33 <= EDIT_QUAL_CUTOFF : continue 
                    else: 
                        edit = (i + offset, ref_nt, sample_nt)
                        if edit in reporter_edit_counts.keys():
                            reporter_edit_counts[edit] += 1
                        else: 
                            reporter_edit_counts[edit] = 1
                gname_to_count[gRNA_name] = reporter_edit_counts

            else:
                # Only count the gRNA counts
                if gRNA_name in gname_to_count.keys(): 
                    gname_to_count[gRNA_name] += 1
                else: 
                    gname_to_count[gRNA_name] = 1

    if not gname_to_count:
        error("No reads assigned to gRNA. check the input.")
    return(gname_to_count)


if __name__ == '__main__':
        REPORTER_LENGTH = 32
        GUIDE_LENGTHS = [19, 20]
        EDIT_QUAL_CUTOFF = 30
        EDIT_START_POS = 1
        EDIT_END_POS = 12
#    try:
        parser = get_input_parser()
        args = parser.parse_args()

        args = check_arguments(args)
        N_READS_INPUT = read_count_match(args.R1, args.R2)

        database_id = get_database_id(args)
        OUTPUT_DIRECTORY='CRISPRessoCount_on_%s' % database_id
        if args.output_folder:
                OUTPUT_DIRECTORY=os.path.join(os.path.abspath(args.output_folder),OUTPUT_DIRECTORY)
        
        _jp=lambda filename: os.path.join(OUTPUT_DIRECTORY,filename) #handy function to put a file in the output directory
        log_filename=_jp('CRISPRessoCount_RUNNING_LOG.txt')
        
        
        try:
            os.makedirs(OUTPUT_DIRECTORY)
            info('Creating Folder %s' % OUTPUT_DIRECTORY)
        except:
            warn('Folder %s already exists.' % OUTPUT_DIRECTORY)
        
        logging.getLogger().addHandler(logging.FileHandler(log_filename))
        with open(log_filename,'w+') as outfile:
            outfile.write('[Command used]:\nCRISPRessoCount %s\n\n[Execution log]:\n' % ' '.join(sys.argv))
        info('Done!')


        # Check read name is aligned to be the same in R1 and R2 files.
        # Optionally filter the reads based on the base quality
        
        filtered_R1_filename = _jp(os.path.basename(args.R1).replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz')
        filtered_R2_filename = _jp(os.path.basename(args.R2).replace('.fastq','').replace('.gz','')+'_filtered.fastq.gz')
        if path.exists(filtered_R1_filename) and path.exists(filtered_R2_filename):
            warn("Using preexisting filtered file")
        else:
            filtered_R1_filename, filtered_R2_filename = check_read_names_and_filter_quality(args.R1,
                args.R2, 
                output_R1_filename= filtered_R1_filename,
                output_R2_filename= filtered_R2_filename,
                min_bp_quality=args.min_average_read_quality,
                min_single_bp_quality=args.min_single_bp_quality, 
                qend_R1 = args.qend_R1, 
                qend_R2 = args.qend_R2)

        info('Done!')
        
        N_READS_AFTER_PREPROCESSING = read_count_match(filtered_R1_filename, filtered_R2_filename)
        
        if N_READS_AFTER_PREPROCESSING == 0:             
            raise NoReadsAfterQualityFiltering('No reads in input or no reads survived the average or single bp quality filtering.')
        else:
            info('Number of reads in input:%d\tNumber of reads after filtering:%d\n' % (N_READS_INPUT, N_READS_AFTER_PREPROCESSING))
 
        guides_dict = get_sgRNA_dict_simple(args.sgRNA_file)

        if args.count_reporter: 
            guide_to_reporter = get_guide_info(args.sgRNA_file)
            gRNA_count = count_masked_guides(R1_filename = filtered_R1_filename, 
                R2_filename = filtered_R2_filename, 
                guides_dict = guides_dict, 
                guide_start = args.guide_start, 
                guide_bc_len = args.guide_bc_len, 
                write_nomatch = True, 
                count_reporter_edits = True, 
                guide_info_df = guide_to_reporter)
        else:
            gRNA_count = count_masked_guides(R1_filename = filtered_R1_filename, 
                R2_filename = filtered_R2_filename, 
                guides_dict = guides_dict, 
                guide_start = args.guide_start, 
                guide_bc_len = args.guide_bc_len, 
                write_nomatch = True, 
                count_reporter_edits = False, 
                )

        info('Done!')    

        info('Writing output table...')

        if args.count_reporter:
            df_guide_counts = pd.DataFrame.from_records(
                [(name, pos, ref_base, sample_base, count)
                for name, edit_to_count in gRNA_count.items()
                for (pos, ref_base, sample_base), count in edit_to_count.items() 
                ],
                columns = ["name", "pos", "ref_base", "sample_base", "count"]
            )
            df_guide_counts.sort_values(by = "count", ascending = False).to_csv(
                _jp('readCount_edit_{}.txt'.format(database_id)), sep='\t', index = False)

        else:
            df_guide_counts = pd.Series(gRNA_count, name = "read_counts").to_frame()
            df_guide_counts.index.name='guide_name'
            df_guide_counts['read_%']=df_guide_counts['read_counts']/N_READS_AFTER_PREPROCESSING*100
            df_guide_counts['RPM']=df_guide_counts['read_counts']/N_READS_AFTER_PREPROCESSING*1000000
            df_guide_counts.sort_values(by='read_counts',ascending=False).to_csv(
                _jp('readCount_guide_{}.txt'.format(database_id)), sep='\t')
        
        info('Done!')
        
        if not args.keep_intermediate:
            info('Removing Intermediate files...')
            
            files_to_remove=[]
                             
            if args.min_average_read_quality>0 or args.min_single_bp_quality>0:
                files_to_remove.extend([filtered_R1_filename, filtered_R2_filename])
                
            for file_to_remove in files_to_remove:
                 try:
                     if os.path.islink(file_to_remove):
                         os.unlink(file_to_remove)
                     else:                             
                         os.remove(file_to_remove)
                 except:
                     warn('Skipping:%s' %file_to_remove)    

 
        info('All Done!')
        print(r'''
              )                                             )
             (           ________________________          (
            __)__       | __   __            ___ |        __)__
         C\|     \      |/  ` /  \ |  | |\ |  |  |     C\|     \
           \     /      |\__, \__/ \__/ | \|  |  |       \     /
            \___/       |________________________|        \___/
        ''')
        sys.exit(0)
    
#    except Exception as e:
#        error('\n\nERROR: %s' % e)
#        sys.exit(-1)

                if gRNA_name in gname_to_count.keys(): 
                    gname_to_count[gRNA_name] += 1
                else: 
                    gname_to_count[gRNA_name] = 1
