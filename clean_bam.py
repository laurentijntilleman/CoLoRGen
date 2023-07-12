#!/usr/bin/env python3

# packages
import logging
import argparse
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if __name__ == "__main__":
    # input
    parser = argparse.ArgumentParser(description='get second alignments')
    parser.add_argument('--log_level', metavar='log_level', dest='log_level',
                        type=str, help='log level')
    parser.add_argument('--log_file', metavar='log_file', dest='log_file',
                        type=str, help='log file')
    parser.add_argument('-i','--input', metavar='bam', dest='sample_in',
                        type=str, help='path to the mapped reads (*.bam)')
    parser.add_argument('-o','--output', metavar='prefix', dest='sample_out',
                        type=str, help='path to output the cleand reads (prefix*.bam)')
    args = parser.parse_args()

    # logging
    log_numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(log_numeric_level, int):
        raise ValueError('Invalid log level: %s' % args.log_level)
    logging.basicConfig(level=log_numeric_level, filename=args.log_file,
                        format='%(asctime)s %(message)s')

    logging.info(f'Start clean {args.sample_in}')

    bamfile = pysam.AlignmentFile(args.sample_in)
    out_file = pysam.AlignmentFile(f'{args.sample_out}.bam', "wb", template=bamfile)
    sequences1 = []
    sequences2 = []
    for read in bamfile.fetch():
        if read.mapq != 0:
            out_file.write(read)
            if read.reference_name == 'H1':
                sequences1.append(SeqRecord(seq=Seq(read.seq),
                                            id = read.query_name,
                                            description = "",
                                            letter_annotations = {'phred_quality':list(read.query_qualities)}
                                           ))
            if read.reference_name == 'H2':
                sequences2.append(SeqRecord(seq=Seq(read.seq),
                                            id = read.query_name,
                                            description = "",
                                            letter_annotations = {'phred_quality':list(read.query_qualities)}
                                           ))
    out_file.close()
    SeqIO.write(sequences1,f'{args.sample_out}_H1.fastq','fastq')
    SeqIO.write(sequences2,f'{args.sample_out}_H2.fastq','fastq')
