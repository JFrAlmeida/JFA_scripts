import pandas as pd
import os
import shutil

"""
Takes in a folder with genomes, all with the same fasta format and a csv or similar file with their names and a 
grouping column named bad_group

file columns:
    genome names = "Names"
    group column = "Group"

from it, generates the QUERY_LIST and REFERENCE_LIST necessary for fastany many to many operation

Query list is marked "Yes" in the "Groups" column
Reference list is marked "No" in the "Groups" column
"""

###Code below
input_genomes_folder = "all_together_now"
input_group = "Groups.csv"
output_folder = "fastani_lists"
Query_output_name = "query.txt"
ref_output_name = "ref.txt"


#Make output if it doesn't exist
if os.path.exists(output_folder):
    shutil.rmtree(output_folder)
    os.makedirs(output_folder)
else:
    os.makedirs(output_folder) #if folder dont exist, make it



df_group = pd.read_csv(input_group)

group_yes_df = df_group.loc[df_group["Groups"] == "Yes"]
group_no_df = df_group.loc[df_group["Groups"] == "No"]

query_list = [name_query for name_query in group_yes_df.Names.to_list()]
reference_list = [name_ref for name_ref in group_no_df.Names.to_list()]

export_query = [item + ".fasta" + "\n" for item in query_list]
export_ref = [item + ".fasta" + "\n" for item in reference_list]

#export query
with open(os.path.join(output_folder,Query_output_name), "w") as file:
    for item in export_query:
        file.write(item)
#export reference
with open(os.path.join(output_folder,ref_output_name), "w") as file:
    for item in export_ref:
        file.write(item)
