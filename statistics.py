# %%
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import os 
import pandas
import math

root_fasta_path = "/mnt/c/works/RKI/AMPsConverter/AMP_DB/fasta"
root_img_path = "/mnt/c/works/RKI/AMPsConverter/AMP_DB/img"

def plotCharDistribution(dict_counts, char_totalsize, file_path):
    head, tail = os.path.split(file_path) 
    save_file_path = os.path.join(root_img_path, tail)

    fig = plt.figure()
    fig.suptitle(tail, fontsize=20)
    ax = fig.add_subplot(111)
    # The bar chart, with letters along the horizontal axis and the calculated
    # letter frequencies as percentages as the bar height
    x = range(len(dict_counts))
    letters =  dict_counts.keys()
    ax.bar(x, [val / char_totalsize * 100 for key,val in dict_counts.items()], width=0.8,
        color='g', alpha=0.5, align='center')
    ax.set_xticks(x)
    ax.set_xticklabels(letters)
    ax.tick_params(axis='x', direction='out')
    ax.set_xlim(-0.5, 25.5)
    ax.yaxis.grid(True)
    ax.set_ylabel('Letter frequency, %')
    plt.savefig(save_file_path+'.png')

def countProtien(file):
    # count Sequence 
    records = list(SeqIO.parse(file, "fasta"))
    print("Count Records:", len(records))
    # countAminoAcids
    with open(file) as f:
        counts = Counter()
        for line in f: # iterate over file object, no need to read all contents onto memory
            if line.startswith(">") or line == "": # skip lines that start with > OR  The empty string is a False value.
                continue
            counts.update(line.rstrip())
        total = float(sum(counts.values()))

        for key,val in counts.items():
            print("{}: {}, ({:.2%})".format(key,val, val / total))

    plotCharDistribution(counts, total, file)

# %%
for f in os.listdir(root_fasta_path):
    print( "---------"+ f + "---------")
    file = os.path.join(root_fasta_path, f)
    countProtien(file)
    print( "--------- End of "+ f + "---------")

# %%
# SortedFile
file = "/mnt/c/works/RKI/AMPsConverter/AMP_DB/one_fasta_file.sorted.fasta"
countProtien(file)

# %%
# PlotCluster
file = "/mnt/c/works/RKI/AMPsConverter/AMP_DB/stats/nr100"
countProtien(file)

# %%
# lenghtDistribution
def lengthDistribution(file_path):
    head, tail = os.path.split(file_path) 
    lengths = map(len, SeqIO.parse(file_path, 'fasta'))
    fig = plt.figure(figsize = (50,50))
    ax = fig.gca()
    pandas.Series(lengths).hist(color='blue', bins=200, ax=ax)
    plt.ylabel('Num of Seq with this length ')
    plt.xlabel('Length of Seq')
    plt.yticks(np.arange(0, 10000, 50))
    plt.xticks(np.arange(0, 2500, 50))
    plt.title(tail)
    save_file_path = os.path.join(root_img_path, tail)
    plt.savefig(save_file_path+'.lenDis.png')

file = "/mnt/c/works/RKI/AMPsConverter/AMP_DB/stats/nr100"
lengthDistribution(file)
# %%
# lenghtDistribution For github
def lengthDistribution(file_path):
    head, tail = os.path.split(file_path) 
    lengths = map(len, SeqIO.parse(file_path, 'fasta'))
    fig = plt.figure()
    ax = fig.gca()
    pandas.Series(lengths).hist(color='blue', bins=200, ax=ax)
    plt.ylabel('Num of Seq with this length ')
    plt.xlabel('Length of Seq')
    plt.title(tail)
    save_file_path = os.path.join(root_img_path, tail)
    plt.savefig(save_file_path+'.github.png')

file = "/mnt/c/works/RKI/AMPsConverter/AMP_DB/stats/nr100"
lengthDistribution(file)
# %%
