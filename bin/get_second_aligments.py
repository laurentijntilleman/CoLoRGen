#!/usr/bin/env python3

# packages
import logging
import argparse
import gzip
import pysam
from functions import reverce_complement
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from array import array
from collections import defaultdict

def find_second_aligments(sample_in,sample_out,target_start,target_stop):
    bamfile = pysam.AlignmentFile(sample_in)
    out_file = pysam.AlignmentFile(f'{sample_out}_sec.bam', "wb", template=bamfile)
    split_fastq = f'{sample_out}_split.fastq'
    split_fastq2 = f'{sample_out}_split_sec.fastq'
    sequences = []
    sequences2 = []
    sequence_dict = defaultdict(list)
    read_dict = dict()
    out_file2 = pysam.AlignmentFile(f'{sample_out}_no_sec.bam', "wb", template=bamfile)
    logging.info('Start finding secondary aligments')
    for read in bamfile.fetch():
        if 'HP' in dict(read.tags):
            tag = dict(read.tags)['HP']
        else:
            tag = '0'
        if read.seq: #no secondairy allignments
            if read.reference_start < target_stop or read.reference_start + read.reference_length > target_start:
                sequence_dict[tag].append(SeqRecord(Seq(read.seq),
                                                    id=read.query_name,
                                                    name=read.query_name,
                                                    description="",
                                                    letter_annotations = {'phred_quality':list(read.query_qualities)}))
                read_dict[read.query_name] = SeqRecord(Seq(read.seq),
                                                    id=read.query_name,
                                                    name=read.query_name,
                                                    description="",
                                                    letter_annotations = {'phred_quality':list(read.query_qualities)})
            for cigar in [x for x in read.cigar if x[0] ==4 ]:
                if cigar[1] > 100:
                    if cigar == read.cigar[0]:
                        sequences.append(SeqRecord(seq=Seq(read.seq[0:cigar[1]]),
                                       id = '{}:{}:{}-{}_1_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                       name = '{}:{}:{}-{}_1_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                       description = "",
                                       letter_annotations = {'phred_quality':list(read.query_qualities[0:cigar[1]])}
                                      ))
                    if cigar == read.cigar[-1]:
                        sequences.append(SeqRecord(seq=Seq(read.seq[read.query_length-cigar[1]:]),
                                       id = '{}:{}:{}-{}_2_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                       name = '{}:{}:{}-{}_2_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                       description = "",
                                       letter_annotations = {'phred_quality':list(read.query_qualities[read.query_length-cigar[1]:])}
                                      ))
            out_file2.write(read)
    bamfile.close()
    out_file.close()
    out_file2.close()
    logging.info('Write sequences')
    SeqIO.write(sequences,split_fastq ,'fastq')
    SeqIO.write(sequences2,split_fastq2 ,'fastq')
    for h in sequence_dict:
        SeqIO.write(sequence_dict[h],f'{sample_out}_H{h}.fastq','fastq')
    SeqIO.write([read_dict[x] for x in read_dict],f'{sample_out}.fastq','fastq')


if __name__ == "__main__":
    # input
    parser = argparse.ArgumentParser(description='get second alignments')
    parser.add_argument('--log_level', metavar='log_level', dest='log_level',
                        type=str, help='log level', default= "INFO")
    parser.add_argument('-i','--input', metavar='bam', dest='sample_in',
                        type=str, help='path to the mapped reads')
    parser.add_argument('-o','--output', metavar='prefix', dest='sample_out',
                        type=str, help='output prefix')
    parser.add_argument('-s','--start', metavar='start', dest='target_start',
                        type=int, help='start of the target region')
    parser.add_argument('-t','--stop', metavar='stop', dest='target_stop',
                        type=int, help='stop of the target region')
    args = parser.parse_args()


    # logging
    log_numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(log_numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=log_numeric_level,
                        format='%(asctime)s %(message)s')

    logging.info('Start searching secondary aligments')


    # running program
    find_second_aligments(args.sample_in,args.sample_out,args.target_start,args.target_stop)

    logging.info('End searching secondary aligments')
