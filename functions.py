#!/usr/bin/env python3

# Usful funcitions

def reverce_complement(seq):
    seq_dict = {'A':'T','C':'G','T':'A','G':'C'}
    new_seq = ''.join([seq_dict[x] for x in seq][::-1])
    return new_seq
