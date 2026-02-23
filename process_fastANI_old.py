import pandas   as pd
import os

#to see all columns
pd.set_option('display.max_columns', None)

#Locations
fastani_folder = "/home/jfa/Aquimarina_review/ANIs/fastani/"
fastani_file = fastani_folder + "output_fastANI.csv"
iTOL_order_path = fastani_folder + "iTOL_order.csv"

#Load files
fastani_output = pd.read_csv(fastani_file, sep="\t")
itol_order = pd.read_csv(iTOL_order_path)

#Set up a file containing only the entries of itol_order with ".fasta"
itol_dotfasta = itol_order[itol_order["iTOL order"].str.contains(".fasta")]



#make another dict showing which genomes have .fasta in their name
sorter_dotfasta = dict(zip(itol_dotfasta["iTOL order"], [".fasta"]*len(itol_dotfasta["iTOL order"])))

#add a new column with .fasta for every genome that is suposed to end with .fasta
itol_order["dotfasta"] = itol_order["iTOL order"].apply(lambda x: sorter_dotfasta.get(x)) #now all rows which contain
# ".fasta" in the genome name have a column saying ".fasta"


#Delete ".fasta" from the genome name column for itol order
itol_order.loc[:,"iTOL order"] = itol_order.loc[:,"iTOL order"].str.replace(".fasta","")

#do another dict to slap the .fasta termination in the correct columns of the last fastani table now!!
#im assuming the names are correct now
sorter_dotfasta_right_names = dict(zip(itol_order["iTOL order"], itol_order["dotfasta"]))

#Make a dict with a numerical value for each of the iTOL names in order. Idea is make a column using this dict in the fastani table
itol_order_series = itol_order["iTOL order"]
sorterIndex = dict(zip(itol_order_series, range(len(itol_order_series))))

#calculate alignment factor (mappings over total query fragments
fastani_output["alignment factor"] = fastani_output["bidirectional fragment mappings"]/fastani_output[("total query "
                                                                                                       "fragments")]
#select only comparisons with ANI over 95% and that are not one genome compared to itself
fastani_over95 = (fastani_output.loc[
    (fastani_output["ani"] >= 95) & (fastani_output["Query"] != fastani_output["reference"])
])

#Solve the names in "Query" and "reference" being all fucked up
fastani_over95.loc[:,"Query"] = fastani_over95.loc[:,"Query"].str.replace("/home/microecoevo/JFA/all_together_now/","")
fastani_over95.loc[:,"reference"] = fastani_over95.loc[:,"reference"].str.replace("/home/microecoevo/JFA/all_together_now/","")

#now need to delete all .fasta from the fastani_95
fastani_over95.loc[:,"Query"] = fastani_over95.loc[:,"Query"].str.replace(".fasta","")
fastani_over95.loc[:,"reference"] = fastani_over95.loc[:,"reference"].str.replace(".fasta","")

#Then add column for which the genomes should have .fasta in the name!
fastani_over95.loc[:,"dotfasta"] = fastani_over95.loc[:,"Query"].apply(lambda x: sorter_dotfasta_right_names.get(x))




#populate the fastani table with values from the dict generated from itol_order
fastani_over95.loc[:,"sorter"] = fastani_over95.loc[:,"Query"].apply(lambda x: sorterIndex.get(x))

#add the ".fasta" to the fastani table using the .fasta dict
# sorter_dotfasta_right_names = dict(zip(itol_order["iTOL order"], itol_order["dotfasta"]))
# fastani_over95["dotfasta"] = itol_order["iTOL order"].apply(lambda x: sorter_dotfasta_right_names.get(x)) #now all rows which contain


#Now sort it by ["sorter, "ani"]
fastani_over95_sorted = fastani_over95.sort_values(["sorter", "ani"], ascending=True)

fastani_over95_sorted.to_csv("/home/jfa/Aquimarina_review/ANIs/fastani_over95_sorted.csv", index=False)

#garbage bin
# fastani_output_sortANI = fastani_output.sort_values(by=["ani"], axis=0, ascending=False)
# print(fastani_over95.head())

