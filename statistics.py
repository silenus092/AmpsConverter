# %%
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import os 

root_fasta_path = "/mnt/c/works/RKI/AMPsConverter/AMP_DB/fasta"
root_img_path = "/mnt/c/works/RKI/AMPsConverter/AMP_DB/img"

def plotCharDistribution(dict_counts, char_totalsize, file_path):
    head, tail = os.path.split(file_path) 
    save_file_path = os.path.join(root_img_path, tail)

    fig = plt.figure()
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
            if line.startswith(">") or not line.strip(): # skip lines that start with > OR  The empty string is a False value.
                continue
            counts.update(line.rstrip())
        total = float(sum(counts.values()))
        for key,val in counts.items():
            print("{}: {}, ({:.2%})".format(key,val, val / total))

    plotCharDistribution(counts, total, file)


for f in os.listdir(root_fasta_path):
    print( "---------"+ f + "---------")
    file = os.path.join(root_fasta_path, f)
    countProtien(file)
    print( "--------- End of "+ f + "---------")
# %%
