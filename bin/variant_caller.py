#!/usr/bin/env python3

# packages
import argparse
import pysam
import json
from Bio import SeqIO
import medaka.labels
import medaka.variant
import numpy as np
import pandas as pd
from itertools import groupby
import types
import re
import os

def target_len(cigar_string):
    """
    Given a CIGAR string, return the number of bases consumed from the
    query sequence.
    """
    read_consuming_ops = ("M", "D", "=", "X")
    result = 0
    cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
    for _, length_digits in cig_iter:
        length = int(''.join(length_digits))
        op = next(next(cig_iter)[1])
        if op in read_consuming_ops:
            result += length
    return result


def find_variants(sample_in,con_fasta,output):
    bamfile = pysam.AlignmentFile(sample_in)
    con_fasta_dict = {record.id:record for record in SeqIO.parse(con_fasta,'fasta')}
    haplo_dict = {'H1':{'genes':[],'variants':[]},'H2':{'genes':[],'variants':[]}}
    for read in bamfile.fetch():
        haplo_dict[read.reference_name]['genes'].append(read.qname.split('_')[-1])
        test_seq = read.seq
        ref_seq = str(con_fasta_dict[read.reference_name].seq)
        ref_pos = read.reference_start
        test_pos = 0
        ref_seq_alligned = ''
        test_seq_alligned = ''
        variants = []
        for match in read.cigar:
            if match[0] == 0: #match
                ref_seq_alligned += ref_seq[ref_pos:ref_pos+match[1]]
                test_seq_alligned += test_seq[test_pos:test_pos+match[1]]
                ref_pos += match[1]
                test_pos += match[1]
            elif match[0] == 1: # insertion
                ref_seq_alligned += '-' * match[1]
                test_seq_alligned += test_seq[test_pos:test_pos+match[1]]
                test_pos += match[1]
            elif match[0] == 2: #deletion
                ref_seq_alligned += ref_seq[ref_pos:ref_pos+match[1]]
                test_seq_alligned += '-' * match[1]
                ref_pos += match[1]
        ref_pos = read.reference_start
        test_pos = 0
        convert_table = pd.DataFrame(columns=['con','gene'])
        convert_table_con = []
        convert_table_gene = []
        for i in range(len(ref_seq_alligned)):
            if ref_seq_alligned[i] != test_seq_alligned[i]:
                variants.append([ref_pos,test_pos,ref_seq_alligned[i],test_seq_alligned[i],
                                 ref_seq_alligned[i-2:i+3],test_seq_alligned[i-2:i+3]])
            if ref_seq_alligned[i] != '-':
                ref_pos += 1
            if test_seq_alligned[i] != '-':
                test_pos += 1
            convert_table_con.append(ref_pos)
            convert_table_gene.append(test_pos)
        convert_table['con'] = convert_table_con
        convert_table['gene'] = convert_table_gene
        haplo_dict[read.reference_name]['variants'].append({'gene':read.qname.split('_')[-1],
                                                            'position':read.reference_start,
                                                            'stop':read.reference_start + target_len(read.cigarstring),
                                                            'variants':variants,
                                                            'number_variants':len(variants),
                                                            'convert_table':convert_table})
    all_variants = dict()
    gene_dict = dict()
    for haplotype in haplo_dict:
        gene_dict[haplotype] = []
        all_variants[haplotype] = []
        variant_df = pd.DataFrame()
        num_var = 0
        for variants in haplo_dict[haplotype]['variants']:
            variant_df.loc[0,num_var]=0
            for variant in variants['variants']:
                variant_df.loc[variant[0],num_var]=1
                variant_df.loc[variant[0],f'ref_{num_var}']=variant[3]
                variant_df.loc[variant[0],f'alt_{num_var}']=variant[2]
                variant_df.loc[variant[0],f'pos_{num_var}']=variant[1]
            num_var += 1
        variant_df = variant_df.sort_index().fillna(0)
        gene_start_stop = dict()
        for i in range(num_var):
            gene_start_stop[i]={'min':0,'max':0}
            gene_start_stop[i]['min']=haplo_dict[haplotype]['variants'][i]['position']-1000
            gene_start_stop[i]['max']=haplo_dict[haplotype]['variants'][i]['stop']+1000
        gene_num = -1
        gene_counter = -1
        for position in variant_df.index[1:]:
            genes = []
            variants = np.array([])
            index = variant_df.index.get_loc(position)
            for gene in gene_start_stop:
                if position >= gene_start_stop[gene]['min'] and position <= gene_start_stop[gene]['max']:
                    genes.append(gene)
                    variants = np.append(variants,int(sum(variant_df.iloc[max(index-10,0):index+11][gene])))
            if genes[variants.argmin()] != gene_num:
                if gene_num not in genes:
                    gene_counter += 1
                    gene_dict[haplotype].append(dict(haplo_dict[haplotype]['variants'][genes[variants.argmin()]]))
                elif gene_dict[haplotype][gene_counter]['gene'] != 'hybrid':
                    if sum([position < haplo_dict[haplotype]['variants'][x]['position'] for x in genes]) != 0:
                        gene_dict[haplotype][gene_counter] = dict(haplo_dict[haplotype]['variants'][genes[variants.argmin()]])
                    else:
                        gene_dict[haplotype][gene_counter]['gene0'] = str(gene_dict[haplotype][gene_counter]['gene'])
                        gene_dict[haplotype][gene_counter]['gene1'] = haplo_dict[haplotype]['variants'][genes[variants.argmin()]]['gene']
                        gene_dict[haplotype][gene_counter]['gene'] = 'hybrid'
                        gene_dict[haplotype][gene_counter]['variants0'] = gene_dict[haplotype][gene_counter]['variants']
                        gene_dict[haplotype][gene_counter]['variants1'] = haplo_dict[haplotype]['variants'][genes[variants.argmin()]]['variants']
                        gene_dict[haplotype][gene_counter]['convert_table0'] = gene_dict[haplotype][gene_counter]['convert_table']
                        gene_dict[haplotype][gene_counter]['convert_table1'] = haplo_dict[haplotype]['variants'][genes[variants.argmin()]]['convert_table']
                        gene_dict[haplotype][gene_counter]['breakpoint0'] = haplo_dict[haplotype]['variants'][genes[variants.argmin()]]['position']
                        gene_dict[haplotype][gene_counter]['breakpoint1'] = position
                        gene_dict[haplotype][gene_counter]['stop'] = haplo_dict[haplotype]['variants'][genes[variants.argmin()]]['stop']
                        gene_dict[haplotype][gene_counter]['counter'] = 1
                else:
                    counter = gene_dict[haplotype][gene_counter]['counter'] +1
                    gene_dict[haplotype][gene_counter][f'gene{counter}'] = haplo_dict[haplotype]['variants'][genes[variants.argmin()]]['gene']
                    gene_dict[haplotype][gene_counter][f'variants{counter}'] = haplo_dict[haplotype]['variants'][genes[variants.argmin()]]['variants']
                    gene_dict[haplotype][gene_counter][f'convert_table{counter}'] = haplo_dict[haplotype]['variants'][genes[variants.argmin()]]['convert_table']
                    gene_dict[haplotype][gene_counter][f'breakpoint{counter}'] = position
                    gene_dict[haplotype][gene_counter]['stop'] = haplo_dict[haplotype]['variants'][genes[variants.argmin()]]['stop']
                    gene_dict[haplotype][gene_counter]['counter'] = counter
                gene_num = genes[variants.argmin()]
            variant_df.loc[position,'gene'] = genes[variants.argmin()]
        variant_df.to_csv(f'{output}_{haplotype}.csv')
        for gene in gene_dict[haplotype]:
            if gene['gene'] == 'hybrid':
                for counter in range(gene['counter']):
                    for variant in gene[f'variants{counter}']:
                        if variant[0] >= gene[f'breakpoint{counter}'] and variant[0] < gene[f'breakpoint{counter+1}']:
                            all_variants[haplotype].append(variant)
                for variant in gene[f'variants{counter+1}']:
                    if variant[0] >= gene[f'breakpoint{counter+1}']:
                        all_variants[haplotype].append(variant)
            else:
                all_variants[haplotype].extend(gene['variants'])

    return all_variants,gene_dict


def decode_variants_gene(
        self, sample, ref_seq, gene_variants, allele, ambig_ref=False, return_all=False):
    """Convert network output in sample to a set of `medaka.vcf.Variant` s.

    A consensus sequence is decoded and compared with a reference sequence.
    Both substitution and indel variants that may be multi-base will be
    reported in the output `medaka.vcf.Variant` s.

    :param sample: `medaka.common.Sample`.
    :param ref_seq: reference sequence, should be upper case
    :param ambig_ref: bool, if True, decode variants at ambiguous (N)
        reference positions.
    :param return_all: bool, emit VCF records with `'.'` ALT for reference
        positions which contain no variant.

    :returns: list of `medaka.vcf.Variant` objects.
    """
    logger = medaka.common.get_named_logger('Variants')

    if sample.positions['minor'][0] != 0:
        raise ValueError(
            "The first position of a sample must not be an insertion.")

    # TODO: the code below uses all of: unicode arrays, python strings,
    #  and integer encodings. This could be simplified
    pos = sample.positions
    probs = sample.label_probs
    encoding = self._encoding

    # array of symbols retaining gaps
    predicted = self.decode_consensus(sample, with_gaps=True, dtype='|U1')
    # get reference sequence with insertions marked as '*'
    reference = np.full(len(pos), "*", dtype="|U1")
    reference[pos['minor'] == 0] = np.fromiter(
        ref_seq[pos['major'][0]:pos['major'][-1] + 1],
        dtype='|U1')
    # change reference gene sequence
    ref_pos = pos['major'][0] -1
    t1 = -1
    test_teller = 0
    while test_teller < len(gene_variants) and gene_variants[test_teller][0] < ref_pos:
        test_teller += 1
    other_variants = []
    for nuc in reference:
        if test_teller < len(gene_variants):
            t1 += 1
            if nuc != '*':
                ref_pos += 1
            if gene_variants[test_teller][2] == '-':
                if ref_pos == gene_variants[test_teller][0]-1:
                    if reference[t1+1] == '*':
                        reference[t1+1] = gene_variants[test_teller][3]
                        ref_pos -= 1
                    elif reference[t1+1] == gene_variants[test_teller][3] and reference[t1+2] == '*':
                        reference[t1+2] = gene_variants[test_teller][3]
                        ref_pos -= 1
                    else:
                        logger.info('insert error')
                        other_variants.append([allele]+gene_variants[test_teller])
                        test_teller += 1
                        while test_teller < len(gene_variants) and gene_variants[test_teller][0]-1 == ref_pos:
                            other_variants.append([allele]+gene_variants[test_teller])
                            test_teller += 1
                        test_teller -= 1
                    test_teller += 1
            else:
                if ref_pos == gene_variants[test_teller][0]:
                    if reference[t1] == gene_variants[test_teller][2]:
                        reference[t1] = gene_variants[test_teller][3]
                    else:
                        logger.info('error')
                        other_variants.append([allele]+gene_variants[test_teller])
                    test_teller += 1
    # find variant columns using prescription
    is_variant = self._find_variants(pos['minor'], reference, predicted)

    variants = []
    runs = medaka.common.rle(is_variant)
    runs = runs[np.where(runs['value'])]
    for rlen, rstart, is_var in runs:
        rend = rstart + rlen

        # get the ref/predicted sequence with/without gaps
        var_ref_with_gaps = ''.join(s for s in reference[rstart:rend])
        var_pred_with_gaps = ''.join(s for s in predicted[rstart:rend])
        var_ref = var_ref_with_gaps.replace('*', '')
        var_pred = var_pred_with_gaps.replace('*', '')

        # del followed by insertion can lead to non-variant
        # maybe skip things if reference contains ambiguous symbols
        if var_ref == var_pred and is_var:
            continue
        elif (not ambig_ref and
                not set(var_ref).issubset(set(self.symbols))):
            continue

        # As N is not in encoding, encode N as *
        # This only affects calculation of quals.
        var_ref_encoded = (
            encoding[(s if s != 'N' else '*',)]
            for s in var_ref_with_gaps)
        var_pred_encoded = (
            encoding[(s,)] for s in var_pred_with_gaps)

        # calculate probabilities
        var_probs = probs[rstart:rend]
        ref_probs = np.array(
            [var_probs[i, j] for i, j in enumerate(var_ref_encoded)])
        pred_probs = np.array(
            [var_probs[i, j] for i, j in enumerate(var_pred_encoded)])
        ref_quals = self._phred(1.0 - ref_probs)
        pred_quals = self._phred(1.0 - pred_probs)

        info = dict()
        if self.verbose:
            info = {
                'ref_seq': var_ref_with_gaps,
                'pred_seq': var_pred_with_gaps,
                'ref_qs': ','.join((self._pfmt(q) for q in ref_quals)),
                'pred_qs': ','.join((self._pfmt(q) for q in pred_quals)),
                'ref_q': self._pfmt(sum(ref_quals)),
                'pred_q': self._pfmt(sum(pred_quals)),
                'n_cols': len(pred_quals)}

        # log likelihood ratio
        qual = sum(pred_quals) - sum(ref_quals)
        genotype = {'GT': '1', 'GQ': self._pfmt(qual, 0)}
        # position in reference sequence
        var_pos = pos['major'][rstart]
        if pos['minor'][rstart] != 0:
            # variant starts on insert, prepend ref base (normalization
            # doesnt handle this)
            var_ref = ref_seq[var_pos] + var_ref
            var_pred = ref_seq[var_pos] + var_pred
        # create the variant record and normalize
        variant = medaka.vcf.Variant(
            sample.ref_name, var_pos, var_ref,
            alt=var_pred, filt='PASS', info=info, qual=self._pfmt(qual),
            genotype_data=genotype)
        variant = variant.normalize(reference=ref_seq)
        variants.append(variant)

    if return_all:
        # to avoid complications from deletions and multi-reference spans,
        #  we simply output a record for every minor == 0 position
        sites = pos['minor'] == 0
        _pos = pos['major'][sites]
        _probs = probs[sites]
        _ref = reference[sites]
        _enc = [encoding[(s if s != 'N' else '*',)] for s in _ref]
        _quals = self._phred(
            1.0 - np.array(_probs[np.arange(_probs.shape[0]), _enc]))
        _quals_flt = np.char.mod("%.3f", _quals)
        _quals_int = np.char.mod("%d", np.rint(_quals))
        info = dict()
        for p, base, qf, qi in zip(_pos, _ref, _quals_flt, _quals_int):
            genotype = medaka.vcf.GenotypeData(GT='0', GQ=qi)
            variants.append(
                medaka.vcf.Variant(
                    sample.ref_name, p, base,
                    alt='.', filt='.', info=info, qual=qf,
                    genotype_data=genotype))
        variants.sort(key=lambda x: x.pos)

    return variants,other_variants

class ConvertCoordinates():

    logger = medaka.common.get_named_logger('Variants')
    
    _exon_table = pd.DataFrame()

    def _load_exon_data(self,transcript_name):
        if self._exon_table.empty:
            exon_rows = []
            with open(self.args.ref_gtf) as fh:
                for line in fh:
                    if 'transcript_name' in line and 'exon' in line:
                        line_split = line.split('\t')
                        features = line_split[-1]
                        feature_dict = dict(re.findall(r'(\S+)\s+\"?([\w.-]+)\"?;', features))
                        exon_rows.append([feature_dict['transcript_name'],feature_dict['exon_number'], int(line_split[3]), int(line_split[4])])
            self._exon_table = pd.DataFrame(exon_rows, columns=['transcript_name','exon_number', 'start', 'end'])
        df_gene = self._exon_table[self._exon_table['transcript_name'] == transcript_name]
        exon_dict = dict()
        for row in df_gene.index:
            exon_dict[df_gene.loc[row,'exon_number']] = {'start':df_gene.loc[row,'start'],
                                                         'stop': df_gene.loc[row,'end']}
        return exon_dict

    def _star_alleles(self):
        variants = pd.read_csv(self.args.star_alleles,sep='\t',header=1)
        self.variant_dict = dict()
        for row in variants.index:
            if variants.loc[row,'Variant Start']!= '.':
                variant_start = int(variants.loc[row,'Variant Start'])
                if variant_start not in self.variant_dict:
                    self.variant_dict[variant_start] = {'ref':variants.loc[row,'Reference Allele'],
                                                   'alt':{variants.loc[row,'Variant Allele']:[variants.loc[row,'Haplotype Name']]}}
                elif variants.loc[row,'Variant Allele'] not in self.variant_dict[variant_start]['alt']:
                    self.variant_dict[variant_start]['alt']=variants.loc[row,'Haplotype Name']
                else:
                    self.variant_dict[variant_start]['alt'][variants.loc[row,'Variant Allele']].append(variants.loc[row,'Haplotype Name'])
        self.variant_dict_star = dict()
        for row in variants.index:
            if variants.loc[row,'Variant Start'] != '.':
                variant_haplotype_name = variants.loc[row,'Haplotype Name']
                if variant_haplotype_name not in self.variant_dict_star:
                    self.variant_dict_star[variant_haplotype_name] = [(int(variants.loc[row,'Variant Start']),variants.loc[row,'Reference Allele'],variants.loc[row,'Variant Allele'])]
                else:
                    self.variant_dict_star[variant_haplotype_name].append((int(variants.loc[row,'Variant Start']),variants.loc[row,'Reference Allele'],variants.loc[row,'Variant Allele']))

    def _convert_extra_varaints(self):
        extra_variants = pd.read_csv(f'{self.args.output}_extra.txt',sep='\t')
        position = -1
        new_row = -1
        exta_variants_converted = pd.DataFrame(columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE'])
        for row in extra_variants.index:
            if extra_variants.loc[row,'con_pos'] == position:
                exta_variants_converted.loc[new_row,'REF'] += extra_variants.loc[row,'gene_env'][2].replace('-','')
                exta_variants_converted.loc[new_row,'ALT'] += extra_variants.loc[row,'con_env'][2].replace('-','')
            elif extra_variants.loc[row,'gene_env'][1] == '-':
                if extra_variants.loc[row,'gene_env'][0] != '-':
                    new_row += 1
                    exta_variants_converted.loc[new_row,'#CHROM'] = extra_variants.loc[row,'allele']
                    exta_variants_converted.loc[new_row,'POS'] = extra_variants.loc[row,'con_pos']
                    position = extra_variants.loc[row,'con_pos']
                    exta_variants_converted.loc[new_row,'REF'] = extra_variants.loc[row,'gene_env'][0:3].replace('-','')
                    exta_variants_converted.loc[new_row,'ALT'] = extra_variants.loc[row,'con_env'][0:3].replace('-','')
                else:
                    print('variant error')
                    os.exit(1)
            elif extra_variants.loc[row,'con_env'][1] == '-':
                if extra_variants.loc[row,'con_env'][0] != '-':
                    new_row += 1
                    exta_variants_converted.loc[new_row,'#CHROM'] = extra_variants.loc[row,'allele']
                    exta_variants_converted.loc[new_row,'POS'] = extra_variants.loc[row,'con_pos']
                    position = extra_variants.loc[row,'con_pos']
                    exta_variants_converted.loc[new_row,'REF'] = extra_variants.loc[row,'gene_env'][0:3].replace('-','')
                    exta_variants_converted.loc[new_row,'ALT'] = extra_variants.loc[row,'con_env'][0:3].replace('-','')
                else:
                    print('variant error')
                    os.exit(1)
            else:
                new_row += 1
                exta_variants_converted.loc[new_row,'#CHROM'] = extra_variants.loc[row,'allele']
                exta_variants_converted.loc[new_row,'POS'] = extra_variants.loc[row,'con_pos']
                position = extra_variants.loc[row,'con_pos']
                exta_variants_converted.loc[new_row,'REF'] = extra_variants.loc[row,'gene_env'][1:3].replace('-','')
                exta_variants_converted.loc[new_row,'ALT'] = extra_variants.loc[row,'con_env'][1:3].replace('-','')
        exta_variants_converted['ID'] = '.'
        exta_variants_converted['QUAL'] = '100'
        exta_variants_converted['FILTER'] = 'PASS'
        exta_variants_converted['INFO'] = '.'
        exta_variants_converted['FORMAT'] = 'GT:GQ'
        exta_variants_converted['SAMPLE'] = '1:100'
        exta_variants_converted.to_csv(f'{self.args.output}_extra.vcf',sep='\t',index=False)
        return exta_variants_converted

    def __init__(self,args):
        header_count = 0
        self.args = args
        hap_vcf = f'{self.args.output}.vcf'
        self.header = []
        with open(hap_vcf) as fh:
            for line in fh.readlines():
                if line.startswith("##"):
                    self.header.append(line)
                    header_count += 1
                else:
                    break
        self.variants = pd.concat([pd.read_csv(hap_vcf,sep='\t',header=header_count),self._convert_extra_varaints()],ignore_index=True).sort_values(by=['#CHROM','POS'])
        self._star_alleles()
        with open(args.gene_json) as fh:
            data = json.loads(fh.read())
            self.genes = data['genes']
            self.hybrids = data['hybrids']
            self.star_allele_genes = data['star_allele_genes']

    def _write_vcf(self,variant,counter,row,true_variants,con_allele,len_true_variants):
        with open(f'{self.args.output}_{con_allele}_gene_{counter}_{self.gene_info[con_allele][counter]["gene"]}.vcf','w') as fh:
            for line in self.header:
                if not line.startswith('##contig'):
                    fh.write(line)
            if self.gene_info[con_allele][counter]["gene"] == 'hybrid':
                fh.write(f'##contig=<ID={self.genes[self.gene_info[con_allele][counter]["gene1"]]["chromosome"]},length={self.genes[self.gene_info[con_allele][counter]["gene1"]]["contig_length"]}>\n')
                fh.write(f'##gene={self.gene_info[con_allele][counter]["gene0"]}\n')
                fh.write(f'##gene={self.gene_info[con_allele][counter]["gene1"]}\n')
            else:
                fh.write(f'##contig=<ID={self.genes[self.gene_info[con_allele][counter]["gene"]]["chromosome"]},length={self.genes[self.gene_info[con_allele][counter]["gene"]]["contig_length"]}>\n')
                fh.write(f'##gene={self.gene_info[con_allele][counter]["gene"]}\n')
            fh.write('\t'.join(list(self.variants.columns))+'\n')
            for line in variant:
                fh.write('\t'.join([str(x) for x in line])+'\n')
        break_point_dict = dict()
        if self.gene_info[con_allele][counter]["gene"] == 'hybrid':
            for break_point in range(self.gene_info[con_allele][counter]['counter']):
                convert_table0 = self.gene_info[con_allele][counter][f'convert_table{break_point}']
                gene_start0 = self.genes[self.gene_info[con_allele][counter][f'gene{break_point}']]['start']
                convert_table1 = self.gene_info[con_allele][counter][f'convert_table{break_point+1}']
                gene_start1 = self.genes[self.gene_info[con_allele][counter][f'gene{break_point+1}']]['start']
                break_point_position = self.gene_info[con_allele][counter][f'breakpoint{break_point+1}']
                break_point_dict[break_point]={'gene0':self.gene_info[con_allele][counter][f'gene{break_point}'],
                                               'gene1':self.gene_info[con_allele][counter][f'gene{break_point+1}'],
                                               'position0':int(convert_table0.loc[convert_table0["con"]==break_point_position,"gene"].iloc[0]) + gene_start0  -1,
                                               'position1':int(convert_table1.loc[convert_table1["con"]==break_point_position,"gene"].iloc[0]) + gene_start1  -1}
        if self.variants.loc[row,"#CHROM"] not in self.gene_dict:
            self.gene_dict[con_allele] = {counter:{'gene':self.gene_info[con_allele][counter]["gene"],
                                                   'variants': variant,
                                                   'true_variants':true_variants,
                                                   'len_true_variants':len_true_variants,
                                                   'breakpoints':break_point_dict}}
        else:
            self.gene_dict[con_allele][counter] = {'gene':self.gene_info[con_allele][counter]["gene"],
                                                   'variants': variant,
                                                   'true_variants':true_variants,
                                                   'len_true_variants':len_true_variants,
                                                   'breakpoints':break_point_dict}


    def convert(self,gene_info):
        self.gene_info = gene_info
        self.gene_dict = dict()
        counter = 0
        variant = []
        true_variants = []
        len_true_variants = 0
        breakpoint_number = 0
        breakpoint_number_new = 0
        write = False
        con_allele = self.variants.loc[0,'#CHROM']
        row_allele = self.variants.loc[0,'#CHROM']
        position_num = -1
        for row in self.variants.index:
            if self.variants.loc[row,'POS'] == position_num:
                variant = variant[:-2]
            position_num = self.variants.loc[row,'POS']
            row_allele = self.variants.loc[row,'#CHROM']
            if con_allele != row_allele:
                if write:
                    self._write_vcf(variant,counter,row,true_variants,con_allele,len_true_variants)
                    variant = []
                    true_variants = []
                    len_true_variants = 0
                    breakpoint_number = 0
                    write = False
                counter = 0
            con_allele = row_allele
            if counter < len(self.gene_info[row_allele]):
                if self.gene_info[row_allele][counter]['stop'] < self.variants.loc[row,'POS']:
                    self._write_vcf(variant,counter,row,true_variants,con_allele,len_true_variants)
                    variant = []
                    true_variants = []
                    len_true_variants = 0
                    breakpoint_number = 0
                    write = False
                    counter +=1
            if counter < len(self.gene_info[row_allele]):
                if self.gene_info[row_allele][counter]['position'] <= self.variants.loc[row,'POS']:
                    position = self.variants.loc[row,'POS']
                    if self.gene_info[row_allele][counter]["gene"] == 'hybrid':
                        breakpoint_boal = False
                        for breakpoint_counter in range(self.gene_info[row_allele][counter]['counter']+1):
                            if self.variants.loc[row,'POS'] >= self.gene_info[row_allele][counter][f"breakpoint{breakpoint_counter}"]:
                                if breakpoint_counter > breakpoint_number:
                                    breakpoint_boal = True
                                breakpoint_number_new = breakpoint_counter
                        if breakpoint_boal:
                            convert_table0 = self.gene_info[row_allele][counter][f'convert_table{breakpoint_number_new-1}']
                            gene_start0 = self.genes[self.gene_info[row_allele][counter][f'gene{breakpoint_number_new-1}']]['start']
                            gene0 = self.gene_info[row_allele][counter][f'gene{breakpoint_number_new-1}']
                            convert_table1 = self.gene_info[row_allele][counter][f'convert_table{breakpoint_number_new}']
                            gene_start1 = self.genes[self.gene_info[row_allele][counter][f'gene{breakpoint_number_new}']]['start']
                            gene1 = self.gene_info[row_allele][counter][f'gene{breakpoint_number_new}']
                            variant.append([f'breakpoint: {int(convert_table0.loc[convert_table0["con"]==position,"gene"].iloc[0]) + gene_start0  -1} - {int(convert_table1.loc[convert_table1["con"]==position,"gene"].iloc[0]) + gene_start1  -1}\t{gene0}-{gene1}'])
                        gene = self.gene_info[row_allele][counter][f'gene{breakpoint_number_new}']
                        convert_table = self.gene_info[row_allele][counter][f'convert_table{breakpoint_number_new}']
                        gene_start = self.genes[self.gene_info[row_allele][counter][f'gene{breakpoint_number_new}']]['start']
                        chromosome = self.genes[self.gene_info[row_allele][counter][f'gene{breakpoint_number_new}']]['chromosome']
                        breakpoint_number = breakpoint_number_new
                    else:
                        gene = self.gene_info[row_allele][counter]["gene"]
                        convert_table = self.gene_info[row_allele][counter]['convert_table']
                        gene_start = self.genes[self.gene_info[row_allele][counter]['gene']]['start']
                        chromosome = self.genes[self.gene_info[row_allele][counter]['gene']]['chromosome']
                    self.variants.loc[row,'POS'] = int(convert_table.loc[convert_table['con']==position,'gene'].iloc[0]) + gene_start  -1
                    if gene in self.star_allele_genes:
                        if len(self.variants.loc[row,'REF']) > len(self.variants.loc[row,'ALT']):
                            true_variant_position = self.variants.loc[row,'POS'] + 1
                            alt_allele = '-'
                        elif len(self.variants.loc[row,'REF']) < len(self.variants.loc[row,'ALT']):
                            true_variant_position = self.variants.loc[row,'POS']
                            alt_allele = self.variants.loc[row,'ALT'][1:]
                        else:
                            true_variant_position = self.variants.loc[row,'POS']
                            alt_allele = self.variants.loc[row,'ALT']
                        if true_variant_position in self.variant_dict:
                            if alt_allele in self.variant_dict[true_variant_position]['alt']:
                                len_true_variants += 1
                                true_variants.extend(self.variant_dict[true_variant_position]['alt'][alt_allele])
                            else:
                                self.logger.info(f'allele {alt_allele} do not exists ')
                    self.variants.loc[row,'#CHROM'] = chromosome
                    variant.append(self.variants.loc[row])
                    write = True
            else:
                if write:
                    self._write_vcf(variant,counter,row,true_variants,con_allele,len_true_variants)
                    variant = []
                    true_variants = []
                    len_true_variants = 0
                    breakpoint_number = 0
                    write = False
        if write:
            self._write_vcf(variant,counter,row,true_variants,con_allele,len_true_variants)

    def annotate(self):
        with open(f'{self.args.output}_haplotypes.txt','w') as fh:
            for haplotype in self.gene_dict:
                fh.write(f'{haplotype}\n')
                for gene in self.gene_dict[haplotype]:
                    fh.write(f'Gene {gene} is {self.gene_dict[haplotype][gene]["gene"]} with {len(self.gene_dict[haplotype][gene]["variants"])} variants\n')
                    self.logger.info(self.gene_dict[haplotype][gene]["gene"])
                    if self.gene_dict[haplotype][gene]["gene"] == 'hybrid':
                        gene_intron_exon_dict = dict()
                        start_position = dict()
                        for gene_name in self.genes.keys():
                            start_position[gene_name] = -1
                            gene_intron_exon_dict[gene_name] = {'exons':set(),'introns':set()}
                        for break_point in self.gene_dict[haplotype][gene]["breakpoints"]:
                            break_point_dict = self.gene_dict[haplotype][gene]["breakpoints"][break_point]
                            exon_data = self._load_exon_data(self.genes[break_point_dict['gene0']]['transcript_name'])
                            self.logger.info(break_point_dict)
                            for exon in exon_data:
                                if exon_data[exon]['start'] >= start_position[break_point_dict['gene0']] and exon_data[exon]['start'] <= break_point_dict['position0']:
                                    if exon_data[exon]['stop'] >= start_position[break_point_dict['gene0']] and exon_data[exon]['stop'] <= break_point_dict['position0']:
                                        gene_intron_exon_dict[break_point_dict['gene0']]['exons'].add(int(exon))
                                if str(int(exon)+1) in exon_data:
                                    if exon_data[str(int(exon)+1)]['stop'] >= start_position[break_point_dict['gene0']] and exon_data[str(int(exon)+1)]['stop'] <= break_point_dict['position0']:
                                        if exon_data[exon]['start'] >= start_position[break_point_dict['gene0']] and exon_data[exon]['start'] <= break_point_dict['position0']:
                                            gene_intron_exon_dict[break_point_dict['gene0']]['introns'].add(int(exon)+1)
                                elif exon_data[exon]['start'] >= start_position[break_point_dict['gene0']] and exon_data[exon]['start'] <= break_point_dict['position0']:
                                    gene_intron_exon_dict[break_point_dict['gene0']]['introns'].add(int(exon))
                            start_position[break_point_dict['gene1']] = break_point_dict['position1']
                        exon_data = self._load_exon_data(self.genes[break_point_dict['gene1']]['transcript_name'])
                        for exon in exon_data:
                            if exon_data[exon]['start'] >= start_position[break_point_dict['gene1']]:
                                gene_intron_exon_dict[break_point_dict['gene1']]['exons'].add(int(exon))
                            if exon_data[exon]['stop'] >= start_position[break_point_dict['gene1']]:
                                gene_intron_exon_dict[break_point_dict['gene1']]['exons'].add(int(exon))
                            if str(int(exon)+1) in exon_data:
                                if exon_data[str(int(exon)+1)]['stop'] >= start_position[break_point_dict['gene1']]:
                                    gene_intron_exon_dict[break_point_dict['gene1']]['introns'].add(int(exon)+1)
                            if exon_data[exon]['start'] >= start_position[break_point_dict['gene1']]:
                                gene_intron_exon_dict[break_point_dict['gene1']]['introns'].add(int(exon))
                        for hybrid in self.hybrids:
                            boal_list = []
                            for gene_name in self.genes.keys():
                                for intron_exon in ['introns','exons']:
                                    boal_list.append(set(self.hybrids[hybrid][gene_name][intron_exon]).issubset(gene_intron_exon_dict[gene_name][intron_exon]))
                            if sum(boal_list) == 4:
                                fh.write(f'This hybrid is {hybrid}\n')
                        self.logger.info(gene_intron_exon_dict)
                    else:
                        if self.gene_dict[haplotype][gene]['gene'] in self.star_allele_genes:
                            true_variants = self.gene_dict[haplotype][gene]['true_variants']
                            true_variants_counter = self.gene_dict[haplotype][gene]['len_true_variants']
                            df=pd.DataFrame({'Number': true_variants})
                            df1 = pd.DataFrame(df['Number'].value_counts())
                            df1['Count']=df1['count'].index
                            fh.write(f"{df1['count'].max()} of the {true_variants_counter} star allele variants are of these haplotypes\n")
                            haplotypes = list(df1[df1['count']==df1['count'].max()]['Count'])
                            for haplo in haplotypes:
                                if len(self.variant_dict_star[haplo]) > df1['count'].max():
                                    fh.write(f'{haplo} has more variants than in sample: {len(self.variant_dict_star[haplo])}\n')
                                elif len(self.variant_dict_star[haplo]) == df1['count'].max():
                                    fh.write(f'{haplo} is the correct haplotype\n')
                                else:
                                    fh.write(f'{haplo} has less variants than in sample: {len(self.variant_dict_star[haplo])}\n')



def variants_from_hdf(args):
    """Entry point for variant calling from HDF5 files.

    A `LabelScheme` read from HDF must define both a `decode_variants`
    and `decode_consnesus` method. The latter is used with `join_samples`
    to detect multi-locus variants spanning `Sample` slice boundaries.

    """
    logger = medaka.common.get_named_logger('Variants')

    index = medaka.datastore.DataIndex(args.inputs)

    args.regions = index.regions

    # lookup LabelScheme stored in HDF5
    try:
        label_scheme = index.metadata['label_scheme']
    except KeyError:
        logger.debug(
            "Could not find `label_scheme` metadata in input file, "
            "assuming HaploidLabelScheme.")
        label_scheme = medaka.labels.HaploidLabelScheme()

    label_scheme.decode_variants_gene = types.MethodType(decode_variants_gene, label_scheme)

    logger.debug("Label decoding is:\n{}".format(
        '\n'.join('{}: {}'.format(k, v)
                  for k, v in label_scheme._decoding.items())))

    if not hasattr(label_scheme, 'decode_variants'):
        raise AttributeError(
            '{} does not support decoding of variants'.format(label_scheme))

    if not hasattr(label_scheme, 'decode_consensus'):
        raise AttributeError(
            '{} does not support consensus decoding required '
            'for variant calling.'.format(label_scheme))

    # tell label_scheme whether we want verbose info fields
    label_scheme.verbose = args.verbose

    meta_info = label_scheme.variant_metainfo

    with pysam.FastaFile(args.ref_fasta) as fa:
        lengths = dict(zip(fa.references, fa.lengths))

    all_variants, gene_info = find_variants(args.sample_in,args.ref_fasta,args.output)

    other_variants = []
    with medaka.vcf.VCFWriter(
            f'{args.output}.vcf', 'w', version='4.1',
            contigs=['{},length={}'.format(r.ref_name, lengths[r.ref_name])
                     for r in args.regions],
            meta_info=meta_info) as vcf_writer:
        for reg in args.regions:
            gene_variants = all_variants[reg.ref_name]
            logger.info("Processing {}.".format(reg))
            ref_seq = pysam.FastaFile(args.ref_fasta).fetch(
                reference=reg.ref_name).upper()

            samples = index.yield_from_feature_files([reg])
            trimmed_samples = medaka.common.Sample.trim_samples(samples)
            joined_samples = medaka.variant.join_samples(
                trimmed_samples, ref_seq, label_scheme)

            for sample in joined_samples:
                variants,other_variant = label_scheme.decode_variants_gene(
                    sample, ref_seq, gene_variants, reg.ref_name, ambig_ref=args.ambig_ref,
                    return_all=args.gvcf)
                vcf_writer.write_variants(variants, sort=True)
                other_variants.extend(other_variant)
    with open(f'{args.output}_extra.txt','w') as fh:
        fh.write('allele\tcon_pos\tgene_pos\tcon_allele\tgene_allele\tcon_env\tgene_env')
        for other_variant in other_variants:
            fh.write('\n')
            fh.write('\t'.join([str(x) for x in other_variant]))

    convert_object = ConvertCoordinates(args)
    convert_object.convert(gene_info)
    convert_object.annotate()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('ref_fasta', help='Reference sequence .fasta file.')
    parser.add_argument('sample_in', help='Mapped genes to reference sequence .bam file.')
    parser.add_argument('inputs', help='Consensus .hdf files.')
    parser.add_argument('star_alleles', help='path to the star_allels')
    parser.add_argument('gene_json', help='json file with the hybrids and gene coordinates')
    parser.add_argument('output', help='Output prefix', default='medaka')
    parser.add_argument('ref_gtf', help='Reference gtf', default='ref_gtf')
    parser.add_argument('--verbose', action='store_true',
                            help='Populate VCF info fields.')
    parser.add_argument('--ambig_ref', action='store_true',
                         help='Decode variants at ambiguous reference positions.')
    parser.add_argument('--gvcf', action='store_true',
                         help='Output VCF records for reference loci predicted to be non-variant.')

    args = parser.parse_args()

    variants_from_hdf(args)
