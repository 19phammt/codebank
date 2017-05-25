"""Library related to any kind of file reading/writing, and
some data wrangling tools.
"""

import os
import pandas as pd
import numpy as np

from Bio import SeqIO
# from Bio.SeqFeature import SeqFeature, FeatureLocation


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


def prep_fasta_contents(d0_in, header_cols, seq_col='seq', preview=0):
    """Reads a dataframe into a list of lists, in fasta format.

    Params
    ------
    d0_in: input dataframe
    header_cols: list of str. Whatever columns from d0_in to put in the fasta header
    seq_col: str. The name of the column which contains actual sequence data.
        default = 'seq'.
    preview: int. No. of records to preview after the function is finished.

    Returns
    -------
    contents: a list of lists
    """

    print("WARNING: column selection logic changed: user now specifies columns to place in header!")

    contents = []
    for index, row in d0_in.iterrows():
        header_ls = row[header_cols]
        header_ls = [str(item) for item in header_ls]
        fasta_header = ">" + "|".join(header_ls)
        seq = str(row[seq_col])
        contents.append(fasta_header)
        contents.append(seq)

    # display contents
    if preview > 0:
        for i in range(preview*2):
            if i%2 == 0:# for even rows
                print(contents[i])
            else:
                print(contents[i][:20]+ "..." + contents[i][-20:])

    return contents


def read_flu_data(fn, fmt='fasta'):
    """Reads a file, duh.
    For the moment, this function only works for flu data; that is, the
    name of the virus contains all relevant meta-data, so that we need only
    extract the name and sequence.

    Params
    ------
    fn: str. Path to file.
    fmt: file format. So far, extensions are 'fasta', 'genbank', 'txt'.
    In particular, 'txt' files are still expected to have a fasta-like format.

    Returns
    -------
    seq_df: a pandas dataframe, with colnames:
    Isolate_Id','Isolate_Name','Lineage','Collection_date','seq'
    """
    try:
        handle = open(fn, 'r')
    except IOError:
        print("Can't find the file %s" % fn)
        pass
    else:
        if fmt == 'fasta':
            # This assumes that the input fasta file comes from GISAID, with
            # the header specified as:
            # Isolate ID|Isolate name|Lineage|Collection date
            seq_ls = []
            for record in SeqIO.parse(fn, fmt):
                header_ls = record.id.split("|")
                row = header_ls + [str(record.seq)]
                seq_ls.append(row)

            seq_df = pd.DataFrame(data=seq_ls, columns=
            ['Isolate_Id','Isolate_Name','Lineage','Collection_date','seq'])
        elif fmt == 'txt':
            print('This is code to read .txt files')
        elif fmt == 'genbank':
            print('This is code to read .gb files')


def update_main(input_dict, conn):
    """
    UNUSED. DEPRECATE?

    TO DO: Make the keys of the input dict dynamic, instead of hardcoded.

    DB IO: updates the main table. The columns are static to enforce db discipline.
    By convention, every db (I create) has a 'main' table, where all raw data is.

    Params
    ------
    input_dict: a dictionary with the following keys wth keys = columns, and each
    key only has 1 value associated with it. Altogether, this represents 1 row
    of data to be added to main, so this function will have to be iteratively run over
    many rows of data.
    c = conn.cursor object
    conn = sqlite3.connect() object

    Returns
    -------
    Updates master.db
    """
    c = conn.cursor()
    fieldnames = list(input_dict.keys())
    update_cmd = "INSERT INTO MAIN(" + ', '.join(fieldnames) + ') VALUES (' + ', '.join(['?']*len(fieldnames)) + ')'
    update_vals = []
    for key in fieldnames:
        update_vals.append(input_dict[key])
    update_tuple = tuple(update_vals)

    c.execute(update_cmd, update_vals)
    conn.commit()


def pivot_raw_tbl(df):
    """Pivots a GISAID dataframe of raw data such that each segment is its own column. That is, turns this:

    name_id      | segment | seq
    -------------|---------|----
    roger        |  PB1    | agct...
    roger        |  PB2    | ggtc...
    roger        |  HA     | aatc...

    into this:

    name_id      |   PB1   |   PB2   |   HA
    -------------|---------|---------|--------
    roger        | agct... | ggtc... | aatc...

    The column "name_id" is a concatenation of the isolate name, and its Epiflu ID, delimited by '|'. This is
    necessary as some isolate names can have multiple IDs for some reason (but, thankfully, not the other way
    round.)

    Params
    ------
    df: a pandas dataframe. Must have the 'name_id' and 'segment' columns.

    Returns
    -------
    df_piv: the requested pivoted dataframe.
    """

    #d0a: the partition that will remain unpivoted. d0b: the partition to be pivoted.
    dfa_selection = []; dfb_selection = []
    # select the columns for the partition that will remain unpivoted
    for col in df.columns.values:
        if col not in ['segment', 'seq', 'segment_number']:
            dfa_selection.append(col)

    dfa = df[dfa_selection].drop_duplicates()
    dfb = df[['name_id', 'segment', 'seq']].drop_duplicates()
    dfb_T = dfb.pivot(index='name_id', columns='segment', values='seq')
    dfb_T.columns.name=None
    dfb_T.reset_index(inplace=True)

    df = pd.merge(dfa, dfb_T, on='name_id')
    return df


def get_meta_from_csv(path_meta, d0, meta_cols=["location"], verbose=True):
    """Retrieves metadata from a given csv.
    IMPT: the input df must have name_id as a unique identifier!
    e.g. say d0 =
        name_id      |   PB1   |   PB2   |   HA
        -------------|---------|---------|--------
        sarah        | agct... | ggtc... | aatc...
        elizabeth    | aggt... | ggtc... | aatc...
        jones        | agct... | ggtc... | agtc...

    and the csv =
        name_id      | location |   bday
        -------------|----------|---------
        sarah        | UK/Leeds | ggtc...
        elizabeth    | US/Iowa | ggtc...

    result:
        name_id      |   PB1   |   PB2   |   HA   | location
        -------------|---------|---------|--------|----------
        sarah        | agct... | ggtc... | aatc...| UK/Leeds
        elizabeth    | aggt... | ggtc... | aatc...| US/Iowa
        jones        | agct... | ggtc... | agtc...| *

    PARAMS
    ------
    path_meta: str; path to meta.csv
    d0: input dataframe which needs to have metadata columns added
    meta_cols: columns from meta.csv which you wish to add to d0

    RETURNS
    -------
    d0: the input column, plus the additional requested columns.
    """

    if verbose:
        print("WARNING: only tested for extracting only the 'location' column!")

    try:
        with open(path_meta) as f:
            # Extract metadata from csv
            d0_meta = pd.read_csv(path_meta)
            if verbose:
                print("Meta file located.")

            # Rename certain key columns
            d0_meta.rename(columns={"Isolate_Name":"iso_name",
                                    "Isolate_Id": "iso_id",
                                    "Collection_Date":"cdate",
                                    "Location":"location"},
                           inplace=True)

            # Create a name_id column in d0_meta
            name_id_ls = []
            for index, row in d0_meta.iterrows():
                 name_id_ls.append(row['iso_name'] + '|' + row['iso_id'])
            d0_meta['name_id'] = name_id_ls

            # Select just the columns required, and the key to merge on
            d0_meta = d0_meta[['name_id', 'location']]

            d0 = pd.merge(d0, d0_meta, how='outer', on='name_id')
            f.close()
    except IOError:
        print("%s Metadata file not found!" % path_meta)

    if verbose:
        meta_nid = set(d0_meta['name_id']) # n_uq nids in metadata
        nid = set(d0['name_id']) # n_uq nids in the original fasta
        print("No. of unique name_ids in meta/actual = %s, %s" % (len(meta_nid), len(nid)))
        print("No. of common name_ids = %s " % len(nid.intersection(meta_nid)))
        print("Updated d0.shape = %s " % (d0.shape,))

    return d0


def get_seqs_from_fasta(df):
    """A related operation to get_meta_from_csv(): this time, we add (or override)
    a sequence column to a metadata df, instead of adding metadata columns to
    a seq data df.
    """
    return "WIP"


def location_split(loc_ls, max_loc_len=5, verbose=False):
    """Splits a list of locations, where each location is in the format:
    'continent/country/state/city/district', into an array of 5. Location
    entries with less than 5 elements are padded with a '*' to represent empty
    values.

    Params
    ------
    loc_ls: list of str; a list of locations.
    max_loc_len: int; number of levels of the location hierarchy.

    Returns
    -------
    loc_arr: str array, shape (n_locations, 5).
    """

    loc_arr = []

    counter = 0
    for line in loc_ls:
        ln = line.split("/")
        for i in range(len(ln)):
            ln[i] = ln[i].strip()
        # Pad ln with '*'
        ln = ln + ["*"]*(max_loc_len-len(ln))

        if verbose:
            if counter < 5:
                print(ln)
                counter +=1
        loc_arr.append(ln)

    return np.array(loc_arr)


def adjust_raw_loc(df):
    """Wrapper for location_split(). Formats the existing 'location' column,
    and adds the location array as 5 new columns.
    """
    loc_ls = list(df.location)
    for i in range(len(loc_ls)):
        loca = loc_ls[i].split("/")
        loca = [i.strip() for i in loca]
        loca = "/".join(loca)
        loc_ls[i] = loca
    df["location"] = loc_ls

    loc_arr = location_split(loc_ls)
    df["continent"] = loc_arr[:, 0]
    df["country"] = loc_arr[:, 1]
    df["state"] = loc_arr[:, 2]
    df["city"] = loc_arr[:, 3]
    df["district"] = loc_arr[:, 4]

    return df


def xls_to_csv(fn):
    fn_out = fn[:-4] + ".csv"
    wb = xlrd.open_workbook(fn)
    sheet = wb.sheet_by_name('Tabelle1')
    with open(fn_out, "w") as f:
        wr = csv.writer(f, quoting=csv.QUOTE_ALL)
        for rownum in range(sheet.nrows):
            wr.writerow(sheet.row_values(rownum))


def msa_screen():
    return "WIP"
