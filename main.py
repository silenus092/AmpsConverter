"""
DBAASP Database of Antimicrobial Activity and Structure of Peptides
AVPdb | database of experimentally validated antiviral peptides
HIPdb | A Database of Experimentally Validated HIV Inhibiting Peptide
Hemolytik Database of Hemolytic activity of AMPs.
BaAMPs | Biofilm-active AMPs database: antimicrobial peptides (AMPs) specifically tested against microbial biofilms.
MilkAMP A database for milk AMPs (a comprehensive database of antimicrobial peptides of dairy origin)
DADP: Database of Anuran Defense Peptides
CancerPPD Database of Anticancer peptides and proteins
AntiCP Consists of 225 antimicrobial peptides with anticancer properties.
LAMP2: An update to LAMP database linking antimicrobial peptide.
Antimicrobial Peptide Database (APD3, http://aps.unmc.edu/AP/main.php)
EnzyBase  Database of enzybiotics (lysins, autolysins, lysozymes, and large bacteriocins) (Not sure)
EROP-Moscow The EROP-Moscow oligopeptide database
BACTIBASE : Database Dedicated to Bacteriocin (Filter by antimicrobial keyword  )
DRAMP 2.0, an updated data repository of antimicrobial peptides
AVPpred: collection and prediction of highly effective antiviral peptides. 
PeptideDB database assembles all naturally occurring signalling peptides from animal source
dbAMPv1.4, an integrated resource for exploring antimicrobial peptides
inverpred
Antifp
ADAM
AntiTbPred
"""
# %%
import pandas as pd
import os 

root_path = "/mnt/c/works/RKI/AMPsConverter/AMP_DB/AMPs"
root_output_path = "/mnt/c/works/RKI/AMPsConverter/AMP_DB/fasta"


def createID(ID, source, seq_type=None, description=None):
    header = str(">" + ID + "|" + source)
    if seq_type is not None:
        header = header + "|" + seq_type
    if description is not None:
        header = header + "|" + description
    return header 

# %%
# DRAMP
# 1. remove  "-spacer-"
# 2. remove whitespace

def converterForDRAMP(filename):
    input_path = os.path.join(root_path, filename+".xlsx")
    out_path = os.path.join(root_output_path, filename+".fasta")
    # Load File
    df = pd.read_excel(input_path, header=[0])
    print("size:"+ str(len(df.index)))
    df = groupDRAMPDuplication(df)
    with open(out_path, 'w') as file:
        for index, row in df.iterrows():
            if "unknown" in  row['Sequence'].lower():
                continue
            if  "-spacer-" in  row['Sequence']:
                continue
            header = createID(row['DRAMP_ID'], filename)
            #print(header)
            file.write(header + '\n')
            seq = row['Sequence'] 
            if " " in seq:
                #print("Found whitespace:"+seq)
                seq = seq.replace(" ", "")
            if "-" in seq:
                print("Found '-':"+seq)
            file.write(seq + '\n')

def groupDRAMPDuplication(df):
    df['DRAMP_ID'] = df['DRAMP_ID'].astype(str)
    grouped_df = df.groupby("Sequence")
    grouped_lists = grouped_df["DRAMP_ID"].agg(lambda column: ",".join(column))
    grouped_lists = grouped_lists.reset_index()
    return grouped_lists

file_names=["DRAMP_Antibacterial_amps",
            "DRAMP_Anticancer_amps",
            "DRAMP_Antifungal_amps",
            "DRAMP_Anti-Gram-_amps",
            "DRAMP_Anti-Gram+_amps",
            "DRAMP_Antimicrobial_amps",
            "DRAMP_Antiparasitic_amps",
            "DRAMP_Antiprotozoal_amps",
            "DRAMP_Antitumor_amps",
            "DRAMP_Antiviral_amps",
            "DRAMP_Insecticidal_amps",
            ]

print("--- DRAMP ----")

for file_name in file_names:
    print(file_name)
    converterForDRAMP(file_name)

print("--- End of DRAMP ----")

# %%
# CancerPPD
# lower case?
searchfor = ["not available", "structure given"]

def converterForCancerPPD(filename):
    input_path = os.path.join(root_path, filename+".txt")
    out_path = os.path.join(root_output_path, filename+".fasta")
    # Load File
    df = pd.read_csv(input_path, sep='\t', header=[0])
    df = groupCancerPPDDuplication(df)
    df = df[~df.Sequence.str.contains('|'.join(searchfor), case=False)]

    with open(out_path, 'w') as file:
        for index, row in df.iterrows():
            header = createID(str(row['id']), filename)
            #print(header)
            file.write(header + '\n')
            seq = row['Sequence'].upper() # Convert to upper case
            if " " in seq:
                print("Found whitespace:"+seq)
                seq = seq.replace(" ", "")
            file.write(seq + '\n')

def groupCancerPPDDuplication(df):
    df['id'] = df['id'].astype(str)
    grouped_df = df.groupby("Sequence")
    grouped_lists = grouped_df["id"].agg(lambda column: ",".join(column))
    grouped_lists = grouped_lists.reset_index(name="id")
    return grouped_lists

file_names=["CancerPPD_d_natural",
            #"CancerPPD_d_non-natural",
            "CancerPPD_l_natural",
            #"CancerPPD_l_non-natural",
            "CancerPPD_mix_natural",
            #"CancerPPD_mix_non-natural"
            ]

print("--- CancerPPD ----")

for file_name in file_names:
    print(file_name)
    converterForCancerPPD(file_name)

print("--- End of CancerPPD ----")

# %%
# AVPdb_data
print("--- AVPdb_data ----")
filename = "AVPdb_data"
input_path = os.path.join(root_path, filename+".tsv")
out_path = os.path.join(root_output_path, filename+".fasta")
# Load File
df = pd.read_csv(input_path, header=[0], delimiter='\t')

# Group duplication
df['Id'] = df['Id'].astype(str)
grouped_df = df.groupby("Sequence")
grouped_lists = grouped_df["Id"].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        header = createID(row['Id'], filename)
        # print(header)
        file.write(header + '\n')
        seq = row['Sequence']
        if " " in seq:
            print("Found whitespace:"+seq)
            seq = seq.replace(" ", "")
        # print(seq)
        file.write(seq + '\n')
print("--- End of AVPdb_data ----")

# %%
# BAAMPs_data 
# lower case?
print("--- BAAMPs_data ----")
filename = "BAAMPs_data"
input_path = os.path.join(root_path, filename+".csv")
out_path = os.path.join(root_output_path, filename+".fasta")
# Load File
df = pd.read_csv(input_path, header=[0])
# remove duplicate among same ID
df = df.drop_duplicates(subset=["PeptideID", "PeptideName", "PeptideSequence"], keep="last")
# Group duplication
df['PeptideID'] = df['PeptideID'].astype(str)
grouped_df = df.groupby("PeptideSequence")
grouped_lists = grouped_df["PeptideID"].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(file_name)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        header = createID(row['PeptideID'], filename)
        # print(header)
        file.write(header + '\n')
        seq = row['PeptideSequence'].upper() # Convert to upper case
        # print(seq)
        file.write(seq + '\n')
print("--- End of BAAMPs_data ----")

# %%
# anticp_225_amp_pos
print("--- anticp_225_amp_pos ----")
filename = "anticp_225_amp_pos"
input_path = os.path.join(root_path, filename+".txt")
out_path = os.path.join(root_output_path, filename+".fasta")

# create new ID
df = pd.read_csv(input_path,  names=["Sequence"], sep='\t')
df['Id'] = df.index
df['Id'] = df['Id'].apply(lambda row: str(row)+"_anticp_225")

# Group duplication
df['Id'] = df['Id'].astype(str)
grouped_df = df.groupby("Sequence")
grouped_lists = grouped_df["Id"].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        header = createID(row['Id'], filename)
        # print(header)
        file.write(header + '\n')
        seq = row['Sequence']
        # print(seq)
        file.write(seq + '\n')
print("--- End of anticp_225_amp_pos ----")

# %%
# Hemolytik_allsequences_02_09_2020
# Special char?
# lower case?

print("--- Hemolytik_allsequences ----")
filename = "Hemolytik_allsequences_02_09_2020"
input_path = os.path.join(root_path, filename+".txt")
out_path = os.path.join(root_output_path, filename+".fasta")

df = pd.read_csv(input_path, header=[0], sep='\t', encoding='utf8' )
# create new ID
df['ID'] = df['ID'].apply(lambda row: str(row)+"_Hemolytik_allsequences")

searchfor = ["β", "[", "]","(" , ")", "/", "-", "Ψ","Δ",
 "1","2","3","4","5","6","7","8","9","0"] 
#df = df[~df.SEQ.str.contains('|'.join(searchfor), case=False)]

# Group duplication
df['ID'] = df['ID'].astype(str)
grouped_df = df.groupby("SEQ")
grouped_lists = grouped_df["ID"].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        seq = row['SEQ']
        if any(x in seq for x in searchfor):
            print("special:"+seq)
            continue
        if any(c.islower() for c in seq):
            seq = seq.upper()
            # print("Upper:"+seq)

        if " " in seq:
            print("Found whitespace:"+seq)
            seq = seq.replace(" ", "")
        header = createID(row['ID'], filename)
        # print(header)
        file.write(header + '\n')
        
        # print(seq)
        file.write(seq + '\n')

print("--- End of Hemolytik_allsequences ----")

# %% 
# HIPdb_data
# Format data ? utf-16 -> utf-8 
print("--- HIPdb_data ----")
filename = "HIPdb_data"
input_path = os.path.join(root_path, filename+".csv")
out_path = os.path.join(root_output_path, filename+".fasta")

df = pd.read_csv(input_path, header=0, encoding="utf-8", usecols=['ID', 'SEQUENCE'])

# Group duplication
df['ID'] = df['ID'].astype(str)
grouped_df = df.groupby("SEQUENCE")
grouped_lists = grouped_df["ID"].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        header = createID(row['ID'], filename)
        # print(header)
        file.write(header + '\n')
        seq = row['SEQUENCE']
        if " " in seq:
           print(seq)
           seq = seq.replace(" ", "")

        # print(seq)
        file.write(seq + '\n')
print("--- End of HIPdb_data ----")


# %% 
# dbaasp_peptides
# lower case?
import csv

print("--- dbaasp_peptides ----")
filename = "dbaasp_peptides"
input_path = os.path.join(root_path, filename+".csv")
out_path = os.path.join(root_output_path, filename+".fasta")

df = pd.read_csv(input_path, header=0, usecols=['ID', 'SEQUENCE'])
# create new ID
df['ID'] = df['ID'].apply(lambda row: str(row)+"_dbaasp")

# Group duplication
df['ID'] = df['ID'].astype(str)
grouped_df = df.groupby("SEQUENCE")
grouped_lists = grouped_df["ID"].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        header = createID(row['ID'], filename)
        # print(header)
        file.write(header + '\n')
        seq = row['SEQUENCE'].upper() # Convert to upper case
        if "HBD3" in seq: # temporally fix
           seq = "GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK"
        if " " in seq:
           print("Found whitespace:"+seq)
           seq = seq.replace(" ", "")

        # print(seq)
        file.write(seq + '\n')

print("--- End of dbaasp_peptides ----")


# %%
# DAPD
print("--- DAPD ----")
filename = "DADP"
input_path = os.path.join(root_path, filename+".json")
out_path = os.path.join(root_output_path, filename+".fasta")

df = pd.read_json(input_path)
# split column in to multiple columns by comma 
df = pd.concat([df['aaData'], df['aaData'].apply(pd.Series)],axis=1)
# Group duplication
df[0] = df[0].astype(str)
grouped_df = df.groupby(6)
grouped_lists = grouped_df[0].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()


print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        if row[6] == "/":
            continue
        header = createID(row[0], filename)
        # print(header)
        file.write(header + '\n')
        seq = row[6]
        # print(seq)
        file.write(seq + '\n')

print("--- End of DAPD ----")
# %%
# DAPD
print("--- milkampdb ----")
filename = "milkampdb"
input_path = os.path.join(root_path, filename+".csv")
out_path = os.path.join(root_output_path, filename+".fasta")

df = pd.read_csv(input_path, header=0, usecols=['id', 'Sequence'])
# trim white space 
df['Sequence'] = df['Sequence'].str.strip()
# remove empty value / Not protein symbol
df = df.dropna()
df = df[~df.Sequence.str.contains("Not determined")]
# chop semicolon 
df['Sequence'] = df['Sequence'].str.split(';').str[0]

# Group duplication
grouped_df = df.groupby('Sequence')
grouped_lists = grouped_df['id'].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        seq = row['Sequence']
        if "/" in seq: # temporally fix
           continue
        header = createID(row['id'], filename)
        # print(header)
        file.write(header + '\n')
        # print(seq)
        file.write(seq + '\n')

print("--- End of milkampdb ----")

# %%
# APD3_update2020_release
print("--- APD3_update2020_release ----")
filename = "APD3_update2020_release"
input_path = os.path.join(root_path, filename+".fasta")
out_path = os.path.join(root_output_path, filename+".fasta")

with open(input_path, mode='r') as in_file, \
     open(out_path, mode='w') as out_file:
    for line in in_file:
        if line.startswith(">"): # skip lines that start with > 
            out_file.write(line)
            continue
        removedWhite_line = line.strip() # strip whitespace
        out_file.write(removedWhite_line + '\n')

print("--- End of APD3_update2020_release ----")

# %%
# enzy2
# natural / synthatic
print("--- enzy2 ----")
from Bio import *

filename = "enzy2"
input_path = os.path.join(root_path, filename+".fasta")
out_path = os.path.join(root_output_path, filename+".fasta")

with open(out_path, 'w') as f_out:
    for seq_record in SeqIO.parse(open(input_path, mode='r'), 'fasta'):
        if "-" in  seq_record.seq:
            print('SequenceID = '  + seq_record.id)
            continue
        r=SeqIO.write(seq_record, f_out, 'fasta')
        if r!=1: print('Error while writing sequence:  ' + seq_record.id)

print("--- End of enzy2 ----")


# %%
# LAMP2
# lower case?
print("--- LAMP2 ----")
filename = "LAMP2"
input_path = os.path.join(root_path, filename+".FASTA")
out_path = os.path.join(root_output_path, filename+".fasta")

with open(input_path, mode='r') as in_file, \
     open(out_path, mode='w') as out_file:
    for line in in_file:
        if line.startswith(">"): # skip lines that start with > 
            out_file.write(line.strip()+ '\n')
            continue
        if " " in line:
            line = line.replace(" ", "")
        line = line.upper()
        out_file.write(line)
    out_file.write("\n")

print("--- End of LAMP2 ----")
# %%
# erop
# +: to denote +H2, which is the open N-terminus
# b: for an acetyl residue or other chemical group at the N-terminus
# -: to denote O-, which is the open C-terminus
# z: for an amide bond at the C-terminus
# J: to denote the pyroglutaminyl linkage, formed by an N-terminal glutamine, due to side-chain reaction with the terminal amine residue
# U: for the (occasional) aminoisobutyric acid residue.
#
print("--- erop ----")
filename = "erop"
input_path = os.path.join(root_path, filename+".FASTA")
out_path = os.path.join(root_output_path, filename+".fasta")

with open(input_path, mode='r') as in_file, \
     open(out_path, mode='w') as out_file:
    for line in in_file:
        if line.startswith(">"): # skip lines that start with > 
            out_file.write(line.strip()+ '\n')
            continue

        line = line.replace("+", "")
        line = line.replace("-", "")
        line = line.replace("b", "")
        line = line.replace("z", "")
        if " " in line:
            line = line.replace(" ", "")
        out_file.write(line)

print("--- End of erop ----")


# %% 
# AVPPred 
print("--- AVPPred ----")

filename = "AVPPred"
input_path = os.path.join(root_path, filename+".tsv")
out_path = os.path.join(root_output_path, filename+".fasta")

df = pd.read_csv(input_path, header=0, usecols=['AVP_ID', 'Sequence'], comment='#', sep="\t")
# trim white space 
df['Sequence'] = df['Sequence'].str.strip()
# remove empty value / Not protein symbol
df = df.dropna()

# Group duplication
grouped_df = df.groupby('Sequence')
grouped_lists = grouped_df['AVP_ID'].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        header = createID(row['AVP_ID'], filename)
        # print(header)
        file.write(header + '\n')
        # print(seq)
        file.write(row['Sequence'] + '\n')
print("--- End of AVPPred ----")

# %% PeptideDB
print("--- peptideDB.anti ----")

filename = "peptideDB.anti"

counter=0
input_path = os.path.join(root_path, filename+".FASTA")
out_path = os.path.join(root_output_path, filename+".fasta")
df = pd.DataFrame(columns=['ID','Sequence'])

with open(input_path, mode='r') as in_file:
    for line in in_file:
        if line in ['\n', '\r\n', ""]:
            continue
        ID = str(counter)+"_"+filename
        counter+=1
        new_row = {'ID':ID, 'Sequence':line}
        df = df.append(new_row, ignore_index=True)

# Group duplication
grouped_df = df.groupby('Sequence')
grouped_lists = grouped_df['ID'].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        seq = row['Sequence']

        header = createID(row['ID'], filename)
        # print(header)
        file.write(header + '\n')
        # print(seq)
        file.write(seq )

print("--- End of peptideDB.anti ----")

# %%
# dbAMPv1.4_validated
# lower case?
print("--- dbAMPv1.4_validated ----")
filename = "dbAMPv1.4_validated"
input_path = os.path.join(root_path, filename+".FASTA")
out_path = os.path.join(root_output_path, filename+".fasta")

with open(input_path, mode='r') as in_file, \
     open(out_path, mode='w') as out_file:
    for line in in_file:
        if line.startswith(">"): # skip lines that start with > 
            out_file.write(line.strip()+ '\n')
            continue
        line = line.upper()
        out_file.write(line)

print("--- End of dbAMPv1.4_validated ----")

# %% inverpred
print("--- inverpred  ----")

filename = "inverpred"
input_path = os.path.join(root_path, filename+".txt")
out_path = os.path.join(root_output_path, filename+".fasta")
df = pd.read_csv(input_path, header=0, sep="\t")
# create new ID
df['ID'] = df.index
df['ID'] = df['ID'].apply(lambda row: str(row)+"_"+filename)

# Group duplication
grouped_df = df.groupby('Secuencia')
grouped_lists = grouped_df['ID'].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()
print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        seq = row['Secuencia']
        header = createID(row['ID'], filename)
        # print(header)
        file.write(header + '\n')
        # print(seq)
        file.write(seq + '\n' )

print("--- End of inverpredd ----")


# %%
# Antifp
print("--- Antifp pos_train_ds3  ----")

filename = "pos_train_ds3"
input_path = os.path.join(root_path, filename+".txt")
out_path = os.path.join(root_output_path, filename+".fasta")
df = pd.read_csv(input_path, header=0, sep="\t")
# create new ID
df['ID'] = df.index
df['ID'] = df['ID'].apply(lambda row: str(row)+"_"+filename)

# Group duplication
grouped_df = df.groupby('Sequence')
grouped_lists = grouped_df['ID'].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()
print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        seq = row['Sequence']
        header = createID(row['ID'], filename)
        # print(header)
        file.write(header + '\n')
        # print(seq)
        file.write(seq + '\n' )
print("--- End of Antifp pos_train_ds3 ----")
# %%
# ADAM
print("---  ADAM  ----")

filename = "adam_test"
input_path = os.path.join(root_path, filename+".txt")
out_path = os.path.join(root_output_path, filename+".fasta")
df = pd.read_csv(input_path, header=0, sep="\t")

"""# Group duplication
grouped_df = df.groupby('Sequence')
grouped_lists = grouped_df['ID'].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        seq = row['Sequence']
        header = createID(row['ID'], filename)
        # print(header)
        file.write(header + '\n')
        # print(seq)
        file.write(seq + '\n' )
print("--- End of ADAM ----")
"""

# %%
# AntiTbPred
print("--- antitbpred  ----")

filename = "antitbpred"
input_path = os.path.join(root_path, filename+".txt")
out_path = os.path.join(root_output_path, filename+".fasta")
df = pd.read_csv(input_path, header=0, sep="\t")
# create new ID
df['ID'] = df.index
df['ID'] = df['ID'].apply(lambda row: str(row)+"_"+filename)

# Group duplication
grouped_df = df.groupby('Sequence')
grouped_lists = grouped_df['ID'].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()
print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        seq = row['Sequence']
        header = createID(row['ID'], filename)
        # print(header)
        file.write(header + '\n')
        # print(seq)
        file.write(seq + '\n' )
print("--- End of antitbpred ----")


###########################  Final steps ##########################################
# %%
# Merge all

DIR = '/mnt/c/works/RKI/AMPsConverter/AMP_DB/'
oh = open( os.path.join(DIR,'one_fasta_file.fasta'), 'w')

for f in os.listdir(root_output_path):
    fh = open(os.path.join(root_output_path, f))
    for line in fh:
        oh.write(line)
    fh.close()
oh.close()
print("------ Merging completed ------")

# %%
# Sort in Descending order ( sort them by length, the longest records first)
from Bio import SeqIO
one_fasta_file = '/mnt/c/works/RKI/AMPsConverter/AMP_DB/one_fasta_file.fasta'
output_fasta_file = '/mnt/c/works/RKI/AMPsConverter/AMP_DB/one_fasta_file.sorted.fasta'

records = list(SeqIO.parse(one_fasta_file, "fasta"))
records.sort(key=lambda r: -len(r))
SeqIO.write(records, output_fasta_file, "fasta")
print("------ Sorting completed ------")


# %%
# Run clustering 
from subprocess import Popen, PIPE

output_report = '/mnt/c/works/RKI/AMPsConverter/AMP_DB/stats/nr100'

process = Popen(['/mnt/c/works/RKI/tools/cdhit/cd-hit',
                    '-i', output_fasta_file,
                    '-o', output_report,
                     "-c","1" ], stdout=PIPE, stderr=PIPE)
stdout, stderr = process.communicate()
print("------ stdout ------")
print(stdout)   
print("------ stderr ------")       
print(stderr)                                      

# %%
