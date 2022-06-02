import pandas as pd
import os
import numpy as np
import collections
import shelve

def prep_tvt_from_func_1(validation_frac = 0.05, test_frac = 0.05, 
               bins = ((1,30),(30,500),(500,999999)),
               chunk_len = 150,
               chunk_offset = 50,
               te_seed = 4711 ):
               
    """-------------------------------------------------------------
    prepare training, validation and test data sets
    Read original data from "contigs_func" directory and merge it 
    with labels from "tn.contig.filter.tsv"
    
    Separate the data into training, validation and test data sets. 
    Thereby stratify the data along the length of genome sequences. 
    The size of each set will be determined with the parameters 
    "validation_frac" and "test_frac". "train_frac" will be computed. 
    
    Tokenize the protein families (simply by assigning a number to 
    each unique protein familiy)
    
    Create chunks from this data of size "chunk_len" and with an 
    offset of "chunk_offset". Define these parameters below.
    
    Output:
    
    "df_chunked" (pandas DataFrame)
        contains all chunks with token_ids and labels. 
        The column "set" assigns eauch chunk to test, validation 
        and train sets. 
        Other columns in this dataframe serve for finding the data 
        quickly in the original files
    
    "tokenize_table" (dictionary)
        contains the translation from token to token_ids
            0 = token for padding
            1 - 10 : reserved    
    -------------------------------------------------------------"""
    
    # set parameters
    train_frac = 1 - validation_frac - test_frac
    np.random.seed(te_seed)
    
    # read data into a dictionary
    df = pd.read_csv(os.path.join("..", "data", "original_data", "tn.contig.filter.tsv"), 
                     sep="\t")
    
    col_origin = []
    col_tokens = []
    col_labels = []
    col_length = []
    col_num_TE = []
    col_max_len_TE = []

    for sequence in sorted(os.listdir(os.path.join("..", "data", "original_data", "contigs_func"))):
        df_c = pd.read_csv(os.path.join("..", "data", "original_data", "contigs_func", sequence), 
                           sep = "\t", 
                           names = ("tokens", "strand"), 
                           skiprows = (1) )
        contig_id = sequence.split(".")[0]
        num_TE = 0
        max_len_TE = 0 
        TE_np = np.zeros(df_c.shape[0])
        
        for index, row in df[df["contig_ID"] == int(contig_id)].iterrows():
            num_TE = num_TE + 1
            max_len_TE = max(max_len_TE, row["len"])        
            
            for i in range(row["start"]-1, row["end"]):
                TE_np[i] = 1
        
        tokens = df_c["tokens"].to_list()
        TE = list(TE_np)
        length = len(tokens)
        
        col_origin.append(sequence)
        col_tokens.append(tokens)
        col_labels.append(TE)
        col_length.append(length)
        col_num_TE.append(num_TE)
        col_max_len_TE.append(max_len_TE)
        
    df_a = pd.DataFrame({ "origin": col_origin,
                          "tokens": col_tokens,
                          "labels": col_labels,
                          "length": col_length,
                          "num_TE": col_num_TE,
                          "max_len_TE": col_max_len_TE 
                        })
                        
    # split data into training, validation and test data                    
    df_a["bin"] = "999"
    df_a["set"] = "999"
    
    i = 0
    for bin in bins:
        df_a["bin"].where((df_a["length"] >= bin[0]) & (df_a["length"] < bin[1]), i, inplace = True)
        df_a.loc[df_a["bin"] == i, "set"] = np.random.choice(3,
                                                             df_a["set"][df_a["bin"] == i].count(),  
                                                             p = [train_frac, validation_frac, test_frac])
        i = i + 1
    
    # assign token_ids to tokens
    c = collections.Counter()
    for index, row in df_a.iterrows():
        c.update(row["tokens"])

    tokens = list(c.keys())
    tokenize_table = dict(zip(tokens, list(np.arange(11,len(tokens)+11)))) # reserve 0 to 10 for special token_ids

    df_a["token_ids"] = df_a["tokens"].apply(lambda x: [ tokenize_table[t] for  t in x])
    
    # compute the chunks
    origin = []
    chunk = []
    dset = []
    tokens = []
    token_ids = []
    attention_masks = []
    labels = []

    for index, row in df_a.iterrows():
        l = len(row["token_ids"])
        raw_t = np.zeros(max(l + chunk_offset, chunk_len)).astype(int)
        raw_t[0: l] = row["token_ids"]
        raw_l = np.zeros(max(l + chunk_offset, chunk_len)).astype(int)
        raw_l[0: l] = row["labels"]
        raw_masks = np.zeros(max(l + chunk_offset, chunk_len)).astype(int)
        raw_masks[0:l] = np.ones(l).astype(int)
        raw_tokens = ['pad' for i in range(0, max(l + chunk_offset, chunk_len))]
        raw_tokens[0: l] = row["tokens"]
        chunk_n = 1

        for i in range(0,max(1, l + chunk_offset - chunk_len), chunk_offset):
            origin.append(row["origin"])
            chunk.append(chunk_n)
            dset.append(row['set'])
            token_ids.append(raw_t[i:i+chunk_len])
            labels.append(raw_l[i:i+chunk_len])
            attention_masks.append(raw_masks[i:i+chunk_len])
            tokens.append(raw_tokens[i:i+chunk_len])
            chunk_n = chunk_n + 1


    df_chunked = pd.DataFrame({"origin": origin,
                               "chunk": chunk,
                               "set": dset,
                               "tokens": tokens,
                               "token_ids": token_ids,
                               "attention_masks" : attention_masks,
                               "labels": labels})
    
    df_chunked.loc[df_chunked["set"] == 0, "set"] = "training"
    df_chunked.loc[df_chunked["set"] == 1, "set"] = "validation"
    df_chunked.loc[df_chunked["set"] == 2, "set"] = "test"
    
    return df_chunked, tokenize_table
    
def read_tvt():
    """-------------------------------------------------------------
    read train, validation and test data
    as prepared by prep_trainvalidationtest_1.ipynb

    Output:

    "df_chunked" (pandas DataFrame)
        contains all chunks with token_ids and labels. 
        The column "set" assigns eauch chunk to test, validation and 
        train sets. 
        Other columns in this dataframe serve for finding the data 
        quickly in the original files

    "tokenize_table" (dictionary)
        contains the translation from token to token_ids
            0 = token for padding
            1 - 10 : reserved    
    -------------------------------------------------------------"""
    
    results = shelve.open(os.path.join("..", "data", "prep_trainvalidationtest_1"), 'r')
    df_chunked = results["chunks"]
    tokenize_table = results["tokenize_table"]
    results.close
    
    return df_chunked, tokenize_table