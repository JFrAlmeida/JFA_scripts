import os
import shutil
import time
import pandas as pd
import scipy.stats
from config import path_main_melange



#Hellow!
time_start = time.time()
print(f"starting {os.path.basename(__file__)}...")


#Import all tables
path_index= "/home/jfa/Aquimarina_review/fresh_results/Melange_extra/statistics/tree_order.csv"
path_algae= "/home/jfa/Aquimarina_review/fresh_results/Outputs/Stats_perm/Algae_presence_counter.csv"
path_chitin= "/home/jfa/Aquimarina_review/fresh_results/Outputs/Stats_perm/Chitin_presence_counter.csv"
path_multismash= "/home/jfa/Aquimarina_review/MULTISMASH/output_multismash/good_counts.tsv"

#index ok
df_index = pd.read_csv(path_index,index_col="Tree_order",sep=",")
df_index.columns = ["color"]

#algae ok
df_algae = pd.read_csv(path_algae,index_col=0)

#chitin ok
df_chitin = pd.read_csv(path_chitin,index_col=0)

#multismash
df_multismash = pd.read_csv(path_multismash, sep= "\t", header= 1, index_col= "record")
df_multismash = df_multismash.drop(["description"], axis=1)

def sort_my_index(df):
    index = df.index.to_list()
    index= sorted(index)
    return index

df_list = [df_index,df_algae,df_chitin,df_multismash]
df_list_2 = df_list

for item_1,item_2 in df_list,df_list_2:
    if item_1 != item_2:
        if sort_my_index(item_1) == sort_my_index(item_2):
            print(f"")

"""

Steps:

Get the three tables processed and neat
get another table with the sum of the three categories, chitin, algae, multismash.

check out how the scipy pearson works, though it seemed simple

"""


#Goodbye!
time_finished = time.time()
print(f"{os.path.basename(__file__)} took {time_finished - time_start} seconds to run")
