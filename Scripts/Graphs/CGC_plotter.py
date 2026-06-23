import sys
import os
import time
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt

from dna_features_viewer import GraphicFeature, GraphicRecord

"""
This script produces gene represenations of cazymes gene clusters. takes in one file at a time, feel free to change it 7
into a function and make it iterate over more files if you wish. It has a filter, controlled by targets, and only 
considers gene clusters where at least one gene with one of these annotations exists.

It does not make legends cause I didnt need it to, again, feel free to edit


The basic inputs are: 
    dbcan -> str, path to "cgc_standard_out.tsv", check below
    targets -> lst, list of annotations to filter for
    colors -> dict, general colors for genes example: {"CAZyme": "#ffcccc","TC": "#ccccff","STP": "#ccffcc"}
    color_chitinase -> dict, colors for specific gene annotations, can match targets or be empty. Gets priority 
    over colors

The file the script takes in is "cgc_standard_out.tsv", example:
CGC#	Gene Type	Contig ID	Protein ID	Gene Start	Gene Stop	Gene Strand	Gene Annotation
CGC1	null	CAKLCA010000001.1	CAKLCA010000001.1_303	358863	360032	+	null
CGC1	CAZyme	CAKLCA010000001.1	CAKLCA010000001.1_304	360239	361609	+	CAZyme|GH13_7
CGC2	TC	CAKLCA010000001.1	CAKLCA010000001.1_418	494902	498144	+	TC|1.B.14.6.11
CGC2	CAZyme	CAKLCA010000001.1	CAKLCA010000001.1_419	498283	501522	-	CAZyme|GH74_e1


"""

#Paths to set
dbcan = "/home/jfa/JFA/easy_CGC/aquimarini_Aq135/cgc_standard_out.tsv"

#Set these
targets = ["GH18","GH19","GH20"]
colors = {
    "CAZyme": "#ffcccc",
    "TC": "#ccccff",
    "STP": "#ccffcc",
    "Sulfatase": "#ffffcc",
    "Peptidase": "#8F43FF"
}

color_chitinase = {
    "GH18": "#990000",
    "GH19": "#994C00",
    "GH20": "#CC00CC",
}
#Nomenclature stuff
print('TC: Transporter, TF: Transcription Factor, STP: Signal Transduction Protein')



#Path to .gff file from easy CGC dbcan in active dev

df = pd.read_csv(dbcan,sep="\t").set_index(["CGC#","Gene Annotation"])

pattern = '|'.join(map(re.escape, targets))

result = df.groupby(level='CGC#').filter(
    lambda g: g.index.get_level_values('Gene Annotation')
                      .str.contains(pattern, na=False)
                      .any()
)

for cgc_id, group in result.groupby(level="CGC#"):

    features = []


    for (cgc, annotation), row in group.iterrows():
        if pd.isna(annotation):
            continue

        strand = 1 if row["Gene Strand"] == "+" else -1

        check = 0
        for item in targets:
            if check == 0:
                if item in annotation:
                    color = color_chitinase.get(item,"#DDDDDD")
                    check = 1
                else:
                    color = colors.get(row["Gene Type"], "#DDDDDD")

        feature = GraphicFeature(
            start=int(row["Gene Start"]),
            end=int(row["Gene Stop"]),
            strand=strand,
            color=color,
            label=annotation
        )
        features.append(feature)

    if not features:
        continue

    min_start = min(f.start for f in features)
    max_end = max(f.end for f in features)

    record = GraphicRecord(
        sequence_length=max_end - min_start,
        features=[
            GraphicFeature(
                start=f.start - min_start,
                end=f.end - min_start,
                strand=f.strand,
                color=f.color,
                label=f.label
            )
            for f in features
        ]
    )

    ax, _ = record.plot(
        figure_width=12,
        strand_in_label_threshold=7
    )

    ax.set_title(cgc_id)
    #plt.show()

    ax.figure.savefig(f"{cgc_id}.png", bbox_inches='tight', dpi=300)
