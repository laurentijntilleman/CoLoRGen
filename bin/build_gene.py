#!/usr/bin/env python3

# packages
import argparse
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def build(args):
    with open(args.gene_json) as fh:
        data = json.loads(fh.read())
        target_region = data['target_region']
        genes = data['genes']
    # writing bed and fasta file of the targeted region
    for record in SeqIO.parse(args.fasta_file,'fasta'):
        if record.id == target_region['chromosome']:
            ref_fasta = record
            SeqIO.write(SeqRecord(seq=ref_fasta.seq[target_region['start']-1:target_region['stop']-1],
                                  id='target',
                                  description=f"{target_region['chromosome']}:{target_region['start']}-{target_region['stop']}"),
                        f"{args.output}/target.fasta",'fasta')
    with open(f"{args.output}/target.bed","w") as fh:
        fh.write(f"{target_region['chromosome']}\t{target_region['start']}\t{target_region['stop']}")
    # writing bed and fasta files of the targeted genes
    seqRecords = []
    for gene in genes:
        for record in SeqIO.parse(args.fasta_file,'fasta'):
            if record.id == genes[gene]['chromosome']:
                ref_fasta = record
                seqRecord = SeqRecord(seq=ref_fasta.seq[genes[gene]['start']-1:genes[gene]['stop']-1],
                              id=gene,
                              description=f"{genes[gene]['chromosome']}:{genes[gene]['start']}-{genes[gene]['stop']}")
                SeqIO.write(seqRecord,f"{args.output}/{gene}.fasta",'fasta')
                seqRecords.append(seqRecord)
        with open(f"{args.output}/{gene}.bed","w") as fh:
            fh.write(f"{genes[gene]['chromosome']}\t{genes[gene]['start']}\t{genes[gene]['stop']}")
    SeqIO.write(seqRecords,f"{args.output}/genes.fasta",'fasta')

if __name__ == "__main__":
    # input
    parser = argparse.ArgumentParser(description='make gene sequence file')
    parser.add_argument('-o','--output', metavar='output', dest='output',
                        type=str, help='output file *.fasta')
    parser.add_argument('-r','--reference', metavar='reference', dest='fasta_file',
                        type=str, help='path to the reference file')
    parser.add_argument('-i','--input', metavar='gene_json', dest='gene_json',
                        type=str, help='json file with the hybrids and gene coordinates')
    args = parser.parse_args()

    build(args)
