#! /usr/bin/python3
""" MSA screener
Tired of scrolling endlessly through thousands and thousands of sequences in Genious?
MSA screener prints out a whole bunch of sequence information for you, by reading the
input fasta file!

Usage:
$ python3 msa_screener.py input_fasta.fas 0.8

Don Teng, 24 May 2017
"""

import sys
import pandas as pd
import numpy as np
import argparse

from Bio import SeqIO


""" ============== DEFS ============== """

def replace_chars(char_ls, new_char, my_string):
    """Given my_string, replaces all occurrences of characters in char_ls with
    new_char.
    e.g. my_string = "the_quick?brown//fox!jumps"
    char_ls = ['_', '?', '//', '!']
    replace_chars(char_ls, " ", my_string)
    >> the quick brown fox jumps
    """
    for ch in char_ls:
        my_string = my_string.replace(ch, new_char)

    return my_string


def read_fasta(fn, verbose=False):
    """Reads a fasta into a list of lists, where each list is one record.

    Params
    ------
    fn: string; filename
    verbose: boolean; verbosity.

    Returns
    -------
    contents: a list of lists.
    """
    print_limit = 0
    contents = []
    for record in SeqIO.parse(fn, "fasta"):
        header = record.description
        seq = str(record.seq)
        row = header.split("|") + [seq]
        contents.append(row)
        if verbose:
            if print_limit < 5:
                preview_seq = row[-1][:10] + "..."
                print_row = row[:-1] + [preview_seq]
                print(print_row)
                print_limit +=1

    # strip whitespaces, which GISAID tends to add
    for i in range(len(contents)):
        for j in range(len(contents[i])):
            contents[i][j] = contents[i][j].strip()

    return contents


""" ============== ARGPARSE ============== """

""" ============== PROC ============== """
fn = sys.argv[1]
char_ls = ["-", "n"]
min_len=0.8

# Read fasta into a dataframe
contents = read_fasta(fn)
fasta_cols = []
# Doesn't matter what the fasta header comprises
# only that the last column is named "seq"
for i in range(1,len(contents[0])-1):
    fasta_cols.append("col"+str(i))
fasta_cols = ["iso_name"] + fasta_cols + ["seq"]

df = pd.DataFrame(data=contents, columns=fasta_cols)


# MSA screen
# replace gapped length
seq_ls = list(df["seq"])
seq_len_ls = []
for i in range(len(seq_ls)):
    seq_ls[i] = replace_chars(char_ls, "", seq_ls[i])
    seq_len_ls.append(len(seq_ls[i]))

seq_len_ls = np.array(seq_len_ls)
df["ungapped_seq"] = seq_ls
medi = np.median(seq_len_ls)

print(""" \n============== MSA screener ============== """)
print("no. of sequences = %s" % df.shape[0])
print("Average seq len = %.3f" % np.average(seq_len_ls))
print("Median seq len = %s" % medi)
print("Minimum seq len required = %s" % int(medi*min_len))
print("Min seq len calculated as %s * median seq length" % min_len)
print("")

d_p_contents = []
for index, row in df.iterrows():
    if len(row["ungapped_seq"]) < int(medi*min_len):
        d_p_contents.append([row["iso_name"],
        len(row["ungapped_seq"]),
        len(row["ungapped_seq"])/medi])

d_p = pd.DataFrame(data=d_p_contents, columns=["iso_name",
"len",
"% of median len"])
if d_p.shape[0] > 0:
    print("Short seqs (< min length required)")
    print(d_p)

print(""" \n============== End of output ============== """)
