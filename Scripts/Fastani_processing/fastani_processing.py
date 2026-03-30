import pandas as pd
import os
import shutil

"""
Takes the results of fastani, filters them for the  ANI and AF thresholds specified below
"""

#Edit these two variables
input_file_name = "Galaxy142-[FastANI on dataset 1-141 Output].tsv"
ANI = 99
AF = 0.6

###Code below
if not os.path.exists("Output"):
    os.makedirs("Output") #if folder dont exist, make it

df = pd.read_csv(input_file_name, sep="\t", header=None, names=['query', 'reference', 'ani', 'matches', 'total','af'])
df['af'] = df["matches"]/df["total"]

#eliminate trailing "_fasta_query" and "_fasta_ref"
df['query'] = df['query'].str[:-12]
df['reference'] = df['reference'].str[:-10]

#process the ANI results
df = df[df["query"]!=df["reference"]]
df = df[df["ani"]>ANI]
df = df[df["af"]>AF]

df.to_csv("Output/FastANI_99_processed.csv",index=False)
