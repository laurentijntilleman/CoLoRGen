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

def find_second_aligments(fasta_file,sample_in,sample_out,target_start,target_stop):
    logging.info('Load all sequenced reads')
    fasta = {record.id:record for record in SeqIO.parse(fasta_file,'fastq')}
    logging.info('Reads loaded')
    bamfile = pysam.AlignmentFile(sample_in)
    out_file = pysam.AlignmentFile(f'{sample_out}_sec.bam', "wb", template=bamfile)
    split_fastq = f'{sample_out}_split.fastq'
    split_fastq2 = f'{sample_out}_split_sec.fastq'
    sequences = []
    sequences2 = []
    sequence_dict = defaultdict(list)
    out_file2 = pysam.AlignmentFile(f'{sample_out}_no_sec.bam', "wb", template=bamfile)
    logging.info('Start finding secondary aligments')
    for read in bamfile.fetch():
        if 'HP' in dict(read.tags):
            tag = dict(read.tags)['HP']
        else:
            tag = '0'
        if read.seq: #no secondairy allignments
            if read.reference_start < target_stop or read.reference_start + read.reference_length > target_start:
                sequence_dict[tag].append(fasta[read.query_name])
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
        else:
            new_seq = fasta[read.qname].seq
            query_qualities = array('B',list(fasta[read.qname].letter_annotations['phred_quality']))
            for cigar in [x for x in read.cigar if x[0] ==4 ]:
                if read.is_reverse:
                    new_seq = reverce_complement(fasta[read.qname].seq)
                    query_qualities = array('B',list(fasta[read.qname].letter_annotations['phred_quality'])[::-1])
                    if cigar[1] > 100:
                        if cigar == read.cigar[0]:
                            sequences2.append(SeqRecord(seq=fasta[read.qname].seq[read.query_length-cigar[1]:],
                                           id = '{}:{}:{}-{}_2_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                           name = '{}:{}:{}-{}_2_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                           description = "",
                                           letter_annotations = {'phred_quality':list(fasta[read.qname].letter_annotations['phred_quality'])[read.query_length-cigar[1]:]}
                                          ))
                        if cigar == read.cigar[-1]:
                            sequences2.append(SeqRecord(seq=fasta[read.qname].seq[0:cigar[1]],
                                           id = '{}:{}:{}-{}_1_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                           name = '{}:{}:{}-{}_1_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                           description = "",
                                           letter_annotations = {'phred_quality':list(fasta[read.qname].letter_annotations['phred_quality'])[0:cigar[1]]}
                                          ))

                else:
                    if cigar[1] > 100:
                        if cigar == read.cigar[0]:
                            sequences2.append(SeqRecord(seq=fasta[read.qname].seq[0:cigar[1]],
                                           id = '{}:{}:{}-{}_1_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                           name = '{}:{}:{}-{}_1_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                           description = "",
                                           letter_annotations = {'phred_quality':list(fasta[read.qname].letter_annotations['phred_quality'])[0:cigar[1]]}
                                          ))
                        if cigar == read.cigar[-1]:
                            sequences2.append(SeqRecord(seq=fasta[read.qname].seq[read.query_length-cigar[1]:],
                                           id = '{}:{}:{}-{}_2_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                           name = '{}:{}:{}-{}_2_H{}'.format(read.query_name,read.reference_name,read.reference_start,read.reference_start+read.reference_length,tag),
                                           description = "",
                                           letter_annotations = {'phred_quality':list(fasta[read.qname].letter_annotations['phred_quality'])[read.query_length-cigar[1]:]}
                                          ))
            write_read = pysam.AlignedSegment()
            write_read.query_name = read.qname
            write_read.query_sequence = str(new_seq)
            write_read.flag = 0
            write_read.reference_id = read.reference_id
            write_read.reference_start = read.reference_start
            write_read.mapping_quality = read.mapping_quality
            write_read.cigar = read.cigar
            write_read.next_reference_id = read.next_reference_id
            write_read.next_reference_start = read.next_reference_start
            write_read.template_length = read.template_length
            write_read.query_qualities = query_qualities
            write_read.tags = read.tags
            out_file.write(write_read)
    bamfile.close()
    out_file.close()
    out_file2.close()
    logging.info('Write sequences')
    SeqIO.write(sequences,split_fastq ,'fastq')
    SeqIO.write(sequences2,split_fastq2 ,'fastq')
    for h in sequence_dict:
        SeqIO.write(sequence_dict[h],f'{sample_out}_H{h}.fastq','fastq')


if __name__ == "__main__":
    # input
    parser = argparse.ArgumentParser(description='get second alignments')
    parser.add_argument('--log_level', metavar='log_level', dest='log_level',
                        type=str, help='log level')
    parser.add_argument('--log_file', metavar='log_file', dest='log_file',
                        type=str, help='log file')
    parser.add_argument('-f','--raw_reads', metavar='raw_reads', dest='fasta_file',
                        type=str, help='path to the raw reads')
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
        raise ValueError('Invalid log level: %s' % args.log_level)
    logging.basicConfig(level=log_numeric_level, filename=args.log_file,
                        format='%(asctime)s %(message)s')

    logging.info('Start searching secondary aligments')


    # running program
    if args.fasta_file.endswith('.gz'):
        fasta_file = gzip.open(args.fasta_file,'rt')
    else:
        fasta_file = args.fasta_file
    find_second_aligments(fasta_file,args.sample_in,args.sample_out,args.target_start,args.target_stop)

    logging.info('End searching secondary aligments')
