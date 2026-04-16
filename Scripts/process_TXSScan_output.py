import os
import shutil

import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def ojoin(input1,input2):
    a = os.path.join(input1,input2)
    return a

"Top level variables, change these to your folder names if you want, or the folders to these names:"
input = "./output"
file_wanted = "best_solution.tsv"

join_df = []
for item in os.listdir(input):
    best_solution_path = ojoin(input,ojoin(item, file_wanted))

    df = pd.read_csv(best_solution_path,
                     sep="\t",
                     skiprows=4)
    df =df.dropna(axis=1,how="all")

    #prepare the columns to make a new df
    df["secretory_system"] = df["hit_gene_ref"].str.split("_").str[0]
    df["completeness"] = np.select(
        [
            df["sys_wholeness"] == 1,
         df["sys_wholeness"] < 1
        ],
        [
            "full",
         df["sys_wholeness"]
        ],
        default= "weird, check manually"
    )

    #new dataframe with just genome name, what is the SS and what is its completeness
    summary_df = df[["replicon","secretory_system","completeness"]]
    summary_df =summary_df.drop_duplicates()
    join_df.append(summary_df)

#Concatenates all dfs
final = pd.concat(join_df, ignore_index= True)
final = final[final["completeness"] != "weird, check manually"]

#get counts of SS per genome
wide_final = final.pivot_table(index= "replicon",
                               columns= "secretory_system",
                               aggfunc= "size",
                               fill_value= 0)
#Statistics summary table
summary = wide_final.describe().T[['mean', 'std', 'min', 'max']]
summary['mode'] = wide_final.mode().iloc[0]

#write everything:
final.to_csv("SS_summary_completeness.csv",sep=",",index=False,header=True)
wide_final.to_csv("SS_counts_pergenome.csv",sep=",",index=True,header=True)
summary.to_csv("SS_statistics.csv",sep=",",index=True,header=True)

print(f"script {os.path.basename(__file__)} ran successfully! -- JFA")