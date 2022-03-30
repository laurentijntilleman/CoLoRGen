#!/usr/bin/env python3

# packages
import logging
import argparse
import gzip
import pysam
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_breakpoints(bam_input):
    bamfile = pysam.AlignmentFile(bam_input)
    breakpoint_dict = {'H1':{'del':pd.DataFrame(columns=['count']),'ins':pd.DataFrame(columns=['count'])},
                       'H2':{'del':pd.DataFrame(columns=['count']),'ins':pd.DataFrame(columns=['count'])}}
    for read in bamfile.fetch():
        haplotype = read.query_name.split('_')[-1]
        breakpoint_position = int(read.query_name.split('_')[-2])
        breakpoint_chr = read.query_name.split(':')[1]
        if breakpoint_chr == 'chr22' and read.reference_name == breakpoint_chr:
            if breakpoint_position == 1:
                breakpoint_start = read.reference_start + read.reference_length
                breakpoint_stop = int(read.query_name.split(':')[2].split('-')[0])
            elif breakpoint_position == 2:
                breakpoint_start = int(read.query_name.split(':')[2].split('-')[1].split('_')[0])
                breakpoint_stop = read.reference_start
            else:
                logging.error('Error: no breakpoint!')
                logging.error(read)
            if breakpoint_start < breakpoint_stop:
                break_type = 'del'
            elif breakpoint_stop < breakpoint_start:
                break_type = 'ins'
            else:
                logging.error('Error: no break_type!')
                logging.error(read)
                break
            if haplotype not in breakpoint_dict:
                breakpoint_dict[haplotype]={'del':pd.DataFrame(columns=['count']),'ins':pd.DataFrame(columns=['count'])}
            if f'{breakpoint_chr}:{breakpoint_start}-{breakpoint_stop}' not in breakpoint_dict[haplotype][break_type].index:
                breakpoint_dict[haplotype][break_type].loc[f'{breakpoint_chr}:{breakpoint_start}-{breakpoint_stop}'] = 1
            else:
                breakpoint_dict[haplotype][break_type].loc[f'{breakpoint_chr}:{breakpoint_start}-{breakpoint_stop}'] += 1
    return breakpoint_dict

def build_consensus(breakpoint_dict,fasta_file,target_chr,target_start,target_stop,sample_out):
    logging.info('Load reference genome')
    ref_fasta = {record.id:record for record in SeqIO.parse(fasta_file,'fasta')}
    logging.info('Reads loaded')
    for haplotype in breakpoint_dict:
        logging.debug(f'Haplotype {haplotype}')
        deletions = list(breakpoint_dict[haplotype]['del'].index[breakpoint_dict[haplotype]['del']['count'] > 2])
        if len(deletions) > 1:
            del_ok = []
            for deletion in deletions:
                del_start = deletion.split(':')[1].split('-')[0]
                del_stop = deletion.split('-')[1]
                if del_ok == []:
                    del_ok = [deletion]
                else:
                    for count in range(len(del_ok)):
                        deletion_ok = del_ok[count]
                        del_ok_start = deletion_ok.split(':')[1].split('-')[0]
                        del_ok_stop = deletion_ok.split('-')[1]
                        if del_ok_stop < del_start or del_ok_start > del_stop:
                            del_ok.append(deletion)
                        else:
                            if breakpoint_dict[haplotype]['del'].loc[deletion_ok,'count'] < breakpoint_dict[haplotype]['del'].loc[deletion,'count'] :
                                del_ok[count] = deletion
            logging.debug(del_ok)
        elif len(deletions) == 1:
            del_ok = [deletions[0]]
            logging.debug(del_ok)
        else:
            logging.debug('no deletion')
            del_ok = []
        insertions = list(breakpoint_dict[haplotype]['ins'].index[breakpoint_dict[haplotype]['ins']['count'] > 2])
        if len(insertions) > 1:
            ins_ok = []
            for instertion in insertions:
                ins_start = instertion.split(':')[1].split('-')[0]
                ins_stop = instertion.split('-')[1]
                if ins_ok == []:
                    ins_ok = [instertion]
                else:
                    for count in range(len(ins_ok)):
                        instertion_ok = ins_ok[count]
                        ins_ok_start = instertion_ok.split(':')[1].split('-')[0]
                        ins_ok_stop = instertion_ok.split('-')[1]
                        if ins_ok_stop > ins_start or ins_ok_start < ins_stop:
                            ins_ok.append(instertion)
                        else:
                            if breakpoint_dict[haplotype]['ins'].loc[instertion_ok,'count'] < breakpoint_dict[haplotype]['ins'].loc[instertion,'count'] :
                                ins_ok[count] = instertion
            logging.debug(ins_ok)
        elif len(insertions) == 1:
            ins_ok = [insertions[0]]
            logging.debug(ins_ok)
        else:
            logging.debug('no insertion')
            ins_ok = []
        ok = pd.DataFrame(columns=['chr','start','stop'])
        count = 0
        for deletion in del_ok:
            ok.loc[count,'chr']=deletion.split(':')[0]
            ok.loc[count,'start']=int(deletion.split(':')[1].split('-')[0])
            ok.loc[count,'stop']=int(deletion.split('-')[1])
            count+=1
        for ins in ins_ok:
            ok.loc[count,'chr']=ins.split(':')[0]
            ok.loc[count,'start']=int(ins.split(':')[1].split('-')[0])
            ok.loc[count,'stop']=int(ins.split('-')[1])
        ok.sort_values('start')
        start = target_start
        seq = Seq('')
        for count in ok.index:
            if ok.loc[count,'chr'] == target_chr:
                seq += ref_fasta[target_chr].seq[start:ok.loc[count,'start']]
                start = ok.loc[count,'stop']
        seq += ref_fasta[target_chr].seq[start:target_stop]
        logging.info('Write sequence')
        description = f'Reference sequence {target_chr}:{target_start}-{target_stop}: haplotype {haplotype}'
        if del_ok != []:
            description += f' with {len(del_ok)} deltions'
        if ins_ok != []:
            description += f' with {len(ins_ok)} insertions'
        SeqIO.write(SeqRecord(seq=seq,id=haplotype,description=description),f'{sample_out}_{haplotype}.fasta','fasta')

if __name__ == "__main__":
    # input
    parser = argparse.ArgumentParser(description='get second alignments')
    parser.add_argument('--log_level', metavar='log_level', dest='log_level',
                        type=str, help='log level')
    parser.add_argument('--log_file', metavar='log_file', dest='log_file',
                        type=str, help='log file')
    parser.add_argument('-i','--input', metavar='bam', dest='sample_in',
                        type=str, help='path to the mapped reads bam')
    parser.add_argument('-o','--output', metavar='prefix', dest='sample_out',
                        type=str, help='output prefix')
    parser.add_argument('-r','--reference', metavar='reference', dest='fasta_file',
                        type=str, help='path to the reference file')
    parser.add_argument('-c','--chromosoom', metavar='chr', dest='target_chr',
                        type=str, help='chromosoom of the target region')
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

    logging.info('Start find breakpoints')
    # running program
    breakpoint_dict = get_breakpoints(args.sample_in)
    with open(f'{args.sample_out}_breakpoints.txt','w') as fh:
        fh.write('Break points of the different haplotypes')
        for h in breakpoint_dict:
            fh.write(f'\n\nHaplotype: {h}')
            for break_type in breakpoint_dict[h]:
                fh.write(f'\n\n{break_type}:\n\n')
                fh.write(breakpoint_dict[h][break_type].to_string())
    logging.debug(breakpoint_dict)
    if args.fasta_file.endswith('.gz'):
        fasta_file = gzip.open(args.fasta_file,'rt')
    else:
        fasta_file = args.fasta_file
    build_consensus(breakpoint_dict,fasta_file,args.target_chr,args.target_start,args.target_stop,args.sample_out)
