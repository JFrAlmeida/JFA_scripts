import pandas as pd
import re

if __name__ == "__main__":
    #paths
    path_in = "/home/jfa/fastani/Galaxy142-[FastANI on dataset 1-141 Output].tsv" #where input file?
    path_out = "/home/jfa/fastani/matrix_output.tsv" #where put output file

pd.set_option('display.max_columns', 500)


#each ahs the order by which I want to reorded my matrix
query_path= "order_fastani/order_fasta_query.tsv"
ref_path = "order_fastani/order_fasta_ref.tsv"

df = pd.read_csv(path_in, sep="\t", header=None, names=['q', 'r', 'ani', 'm', 't'])
#Select ANI>95%
df = df[df["ani"]>95]

matrix = df.pivot(index='q', columns='r', values='ani')
print(matrix.index.to_list())

matrix.index = matrix.index.str.replace(r"(?:GCA|GCF)_[0-9]{9}_[0-9]_","",regex= True)
print("matrix_index mod: ",matrix.index)

matrix.columns = matrix.columns.str.replace(r"(?:GCA|GCF)_[0-9]{9}_[0-9]_","",regex= True)
print("columns mod: ",matrix.columns)
matrix.to_csv(path_out, sep='\t')


#
# query_index = pd.read_csv(query_path)
# ref_columns = pd.read_csv(ref_path



