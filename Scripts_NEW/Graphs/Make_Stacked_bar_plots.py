import os
import shutil
import time
import numpy as np
import pandas as pd
import scipy.stats

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

#Hellow!
time_start = time.time()
print(f"starting {os.path.basename(__file__)}...")

path_import = "output_multismash/good_counts.tsv"
tree_order = "/home/jfa/Aquimarina_review/fresh_results/Melange_extra/statistics/tree_order.csv"
export_path = "/home/jfa/Aquimarina_review/MULTISMASH/multismash_stacked.svg"

df = pd.read_csv(path_import,
                 sep= "\t",
                 header= 1,
                 index_col= "record")

df = df.drop(["description","total_count"],
             axis=1)

df_index = pd.read_csv(tree_order, index_col="Tree_order")

"""
input df with index as sample names, order_df with the order of samples and index as samples names, optionally also a 
column with hexadecimal colors to apply to each corresponding label at the bottom of the graph;


Kwargs that matplotlib processes can be passed as kwargs at the end, heres a few I use more often than not:
bars: edgecolor, linewidth
labels: fontsize,fontstyle, ha, rotation
graph export: dpi= int, format= "str", 
custom kwargs: graph_title= "string" , color_start= "hexadecimal", color_end= "hexadecimal", ex_path="path/to/save/file.svg";
color_map= "matplotlib map name"

color_start and color_end are hexadecimal kwa, if not passed defaults to tab10 or tab20 or virides depending on sample 
size; ex_path is the path to save the file as, include the format termination you want the file to have

If not set these default to:
graph_title: "Graph Title"
loc= "upper right"
dpi= 300
format= "svg"
fontsize= 10
fontstyle= "italic"
ha= "right"
rotation= 45
edgecolor= "0.3"
linewidth= 0.7
width= 0.4

"""

def make_stacked_barplot(input_df, export_path, order_df=None, **kwargs):

    if order_df.empty:
        df_reindexed = input_df
        print("no order_df detected, proceeding unordered")
    else:
        if not sorted(input_df.index.to_list()) == sorted(order_df.index.to_list()):
            print(
                "Your input dataframe and order dataframe indexes are not equal, please fix them and run again! \n Quitting...")
            quit()
        df_reindexed = input_df.reindex(order_df.index)

    height_dict = {}
    for column in df_reindexed.columns:
        clist = df_reindexed.loc[:, column]
        cdict = {
            column: clist
        }
        height_dict.update(cdict)

    #Color maps
    color_start = kwargs.get('color_start', None)
    color_end = kwargs.get('color_end', None)
    my_map = kwargs.get('color_map', None)

    cat_len = len(height_dict)
    if color_start and color_end:
        cmap_colors = [color_start,color_end]
        cmap = LinearSegmentedColormap.from_list("custom_gradient", cmap_colors)
        colors = cmap(np.linspace(0, 1, cat_len))
    elif my_map:
        cmap = plt.get_cmap(my_map)
        colors = cmap(np.linspace(0, 1, cat_len))
    else:
        if cat_len <= 10:
            cmap = plt.cm.tab10
            colors = cmap(np.linspace(0, 1, cat_len))
        elif cat_len > 10 and cat_len <= 20:
            cmap = plt.cm.tab20
            colors = cmap(np.linspace(0, 1, cat_len))
        else:
            cmap = plt.cm.viridis
            colors = cmap(np.linspace(0, 1, cat_len))

    # initialise the graph
    fig, ax = plt.subplots(figsize=(30, 6), constrained_layout=True)

    #build graph
    bottom = np.zeros(len(df_reindexed.index))
    eclr = kwargs.get('edgecolor', "0.3")
    lwdth = kwargs.get('linewidth', 0.7)
    stacked_width = kwargs.get('width', 0.4)  # Of the bars

    #Make the stacked bars
    for key in height_dict:
        height_dict[key] = pd.to_numeric(height_dict[key], errors='coerce').fillna(0)

    for (category, count), color in zip(height_dict.items(), colors):
        ax.bar(df_reindexed.index,
               count,
               width=stacked_width,
               bottom=bottom,
               label=category,
               color=color,
               edgecolor=eclr,
               linewidth=lwdth)
        bottom += count

    # Set the tick labels
    fsize = kwargs.get('fontsize', 10)
    fstyle = kwargs.get('fontstyle', "italic")
    horalign = kwargs.get('ha', "right")
    srotation = kwargs.get('rotation', 45)

    labels_x = np.arange(len(df_reindexed.index))  # the label locations
    ax.set_xticks(ticks=labels_x,
                  labels=df_reindexed.index,
                  minor=False,
                  fontsize=fsize,
                  fontstyle=fstyle,
                  ha=horalign,
                  rotation=srotation)

    # set tick colors based on "Color" columns of tree_order.csv
    if "Color" in order_df.columns:
        for tick_label, color in zip(ax.get_xticklabels(), df_index.Color):
            tick_label.set_color(color)
    else:
        print("Color column not found in tree_order.csv, proceeding with uncoloured labels")

    graph_title = kwargs.get('graph_title', "Graph Title")
    legend_loc = kwargs.get('loc', "upper right")
    dpi = kwargs.get('dpi', 300)
    format = kwargs.get('format', "svg")

    ax.set_title(graph_title)
    ax.legend(loc=legend_loc)
    fig.savefig(fname=export_path,
                dpi=dpi,
                format=format)
    print(f"stacked bar plot printed to: {export_path}")
name = "Secondary metabolism Biosythetic gene clusters of Aquimarina"
make_stacked_barplot(df, export_path, order_df= df_index,color_map="terrain",graph_title=name)

#Goodbye!
time_finished = time.time()
print(f"{os.path.basename(__file__)} took {time_finished - time_start} seconds to run")
