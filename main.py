"""
AVPdb database of experimentally validated antiviral peptides
HIPdb A Database of Experimentally Validated HIV Inhibiting Peptide
AntiBP2 : Server for antibacterial peptide prediction
DBAASP V.3 Database of Antimicrobial Activity and Structure of Peptides
DADP: Database of Anuran Defense Peptides
CancerPPD Database of Anticancer peptides and proteins
AntiTbPred | Prediction of antitubercular peptides
AntiCP Consists of 225 antimicrobial peptides with anticancer properties.
Antifp | is an in silico method, which is developed to predict and design antifungal peptides
ADAM | A Database of Anti-Microbial peptides
Inverpred | is a specialized database of AMPs from invertebrates.
PeptideDB database assembles all naturally occurring signalling peptides from animal source
BACTIBASE : Database Dedicated to Bacteriocin (Filter by antimicrobial keyword )
PeptideDB database assembles all naturally occurring signalling peptides from animal source
Antimicrobial Peptide Database (APD3, https://wangapd3.com/downloads.php)
MilkAMP A database for milk AMPs (a comprehensive database of antimicrobial peptides of dairy origin)
BACTIBASE : Database Dedicated to Bacteriocin (Filter by antimicrobial keyword )
AVPpred: collection and prediction of highly effective antiviral peptides.
ADAPTABLE is both a webserver and data-miner of antimicrobial peptides. 
BaAMPs | Biofilm-active AMPs database: antimicrobial peptides (AMPs) specifically tested against microbial biofilms.
dbAMP v1.4, an integrated resource for exploring antimicrobial peptides
DRAMP 2.0, an updated data repository of antimicrobial peptides
Hemolytik Database of Hemolytic activity of AMPs. <-- เอา ID ของ  CP and toxix ออก
EROP-Moscow The EROP-Moscow oligopeptide database
EnzyBase Database of enzybiotics (lysins, autolysins, lysozymes, and large bacteriocins) (อันที่เป็น ช่องว่างในต่ำเเหน่งที่สองของ ID ถือว่าเป็น synthetic)
AmPEP: Sequence-based prediction of antimicrobial 
-----
Not include:
YADAMP: yet another database of antimicrobial peptides  

"""
# %%
import pandas as pd
import os 
from Bio import SeqIO
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
    input_path = os.path.join(root_path+"/DRAMP/",  filename+".xlsx")
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
            if "-" in row['Sequence']:
                print("Found '-':"+row['Sequence'])
                continue

            header = createID(row['DRAMP_ID'], filename)
            #print(header)
            file.write(header + '\n')
            seq = row['Sequence'] 
            if " " in seq:
                #print("Found whitespace:"+seq)
                seq = seq.replace(" ", "")
            file.write(seq + '\n')

def groupDRAMPDuplication(df):
    # Filter 
    Bad_words = ["Unknown","Not found","Non-antibacterial"
                ,"Non-antimicrobial"]
    df = df[~ df['Activity'].isin(Bad_words)]            
    # Remove empyty 
    df = df[df['Sequence'].notna()]

    df['DRAMP_ID'] = df['DRAMP_ID'].astype(str)
    grouped_df = df.groupby("Sequence")
    grouped_lists = grouped_df["DRAMP_ID"].agg(lambda column: ",".join(column))
    grouped_lists = grouped_lists.reset_index()
    return grouped_lists

file_names = ["general_amps","patent_amps","specific_amps",
              "clinical_amps","DRAMP_Antimicrobial_amps"]


#file_names=["DRAMP_Antibacterial_amps",
#            "DRAMP_Anticancer_amps",
##            "DRAMP_Antifungal_amps",
#            "DRAMP_Anti-Gram-_amps",
#            "DRAMP_Anti-Gram+_amps",
#            "DRAMP_Antimicrobial_amps",
#            "DRAMP_Antiparasitic_amps",
#            "DRAMP_Antiprotozoal_amps",
#            "DRAMP_Antitumor_amps",
#            "DRAMP_Antiviral_amps",
#            "DRAMP_Insecticidal_amps",
#            ]

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
    input_path = os.path.join(root_path + "/CancerPPD/", filename+".txt")
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
input_path = os.path.join(root_path, "BaAMPs/"+filename+".csv")
out_path = os.path.join(root_output_path, filename+".fasta")
# Load File
df = pd.read_csv(input_path, header=0)
# remove duplicate among same ID
df = df.drop_duplicates(subset=["PeptideID", "PeptideName", "PeptideSequence"], keep="last")
# remove any record that wasn't validated 
df = df[df['ExperimentStatus'] ==  "Validated"]
# Group duplication
df['PeptideID'] = df['PeptideID'].astype(str)
grouped_df = df.groupby("PeptideSequence")
grouped_lists = grouped_df["PeptideID"].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
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
        if " " in seq:
            print("Found whitespace:"+seq)
            seq = seq.replace(" ", "")
        # print(seq)
        file.write(seq + '\n')
print("--- End of anticp_225_amp_pos ----")

# %%
# Hemolytik_allsequences_02_09_2020
# Special char?
# lower case?

print("--- Hemolytik_allsequences ----")
file_names=["naturalseqD",
            "naturalseqL",
            "naturalseqmix",
            ]
exclude_file_name = ["CPP", "Toxic"]

df = pd.DataFrame()
ex_df = pd.DataFrame()

for ex_file_name in exclude_file_name:
    input_path = os.path.join(root_path +"/HEMOLYTIK/", ex_file_name+".txt")
    df1 = pd.read_csv(input_path, header=[0], sep='\t', encoding='utf8' )
    ex_df = ex_df.append(df1, ignore_index = True)

for file_name in file_names:
    input_path = os.path.join(root_path+"/HEMOLYTIK/", file_name+".txt")
    df1 = pd.read_csv(input_path, header=[0], sep='\t', encoding='utf8' )
    df = df.append(df1, ignore_index = True)

### start ID excluding 
df = df[~ df.ID.isin(ex_df.ID)]

filename = "Hemolytik"
out_path = os.path.join(root_output_path, filename+".fasta")
# create new ID

df['ID'] = df['ID'].apply(lambda row: str(row)+"_"+filename)
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
# Format data? change utf-16 -> utf-8 
print("--- HIPdb_data ----")
filename = "HIPdb_data"
input_path = os.path.join(root_path, filename+".csv")
out_path = os.path.join(root_output_path, filename+".fasta")

df = pd.read_csv(input_path, sep =',', header=0, encoding="utf-8"
                , usecols=['ID', 'SEQUENCE'])

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
filename = "dbaasp.peptides"
input_path = os.path.join(root_path, filename+".csv")
out_path = os.path.join(root_output_path, filename+".fasta")

df = pd.read_csv(input_path, header=0, usecols=['ID', 'SEQUENCE', 'COMPLEXITY'])
# create new ID
df['ID'] = df['ID'].apply(lambda row: str(row)+"_dbaasp")

print("Original df: "+str(len(df)))
## remove multi-peptides and multimer
df = df[df['COMPLEXITY'] == "Monomer"]
print("Filterd df: "+str(len(df)))
# Group duplication
df['ID'] = df['ID'].astype(str)
grouped_df = df.groupby("SEQUENCE")
grouped_lists = grouped_df["ID"].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()
print("Final df: "+str(len(grouped_lists)))
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

df = df.drop(['aaData'], axis=1)

# Group duplication  , select on column 6 (Bioactive sequence)
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
        if " " in seq:
           print(seq)
           seq = seq.replace(" ", "")

        # print(seq)
        file.write(seq + '\n')

print("--- End of DAPD ----")
# %%
# milkampdb
print("--- milkampdb ----")
filename = "milkampdb"
input_path = os.path.join(root_path, filename+".csv")
out_path = os.path.join(root_output_path, filename+".fasta")

df = pd.read_csv(input_path, header=0, usecols=['id','Activity', 'Sequence'])

# Filter "No activity detected"
df = df[df['Activity'] !=  "No activity detected"]
df = df[df['Activity'] !=  "Not determined"]
# trim white space 
df['Sequence'] = df['Sequence'].str.strip()
# remove empty value / Not protein symbol
df = df.dropna()
df = df[~ df.Sequence.str.contains("Not determined", na=False)]
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
        if "/" in seq: # temporally fix (just ignore it )
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

        splited = seq_record.id.split('|')
        if splited[1] =='':
            continue

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
import xml.etree.ElementTree as et 
filename="LAMP2"
input_path = os.path.join(root_path, filename+".XML")
xtree = et.parse(input_path)
print(xtree)

xroot = xtree.getroot()
print(xroot)


# %%
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

def replace_char_at_index(org_str, index, replacement):
    ''' Replace character at index in string org_str with the
    given replacement character.'''
    new_str = org_str
    if index < len(org_str):
        new_str = org_str[0:index] + replacement + org_str[index + 1:]
    return new_str

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
        # J start with 
        if(line.startswith('J')):
            line = replace_char_at_index(line, 0, '')

        # and end with J 
        if(line.endswith('J')):
            line = replace_char_at_index(line, len(line)-1, '')

        # start with U
        if(line.startswith('U')):
            line = replace_char_at_index(line, 0, '')

        # and end with U 
        if(line.endswith('U')):
            line = replace_char_at_index(line, len(line)-1, '')
    
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
input_path = os.path.join(root_path, filename+".fasta")
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
filename = "InverPep"
input_path = os.path.join(root_path, filename+".csv")
out_path = os.path.join(root_output_path, filename+".fasta")
df = pd.read_csv(input_path, header=0, sep=",")
# create new ID
df['ID'] = df.index
df['ID'] = df['ID'].apply(lambda row: str(row)+"_"+filename)

## Filter Origen
df = df[df["Origen"] == "experimentally validated"]

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
## pos_train_ds1 =  1168 antifungal peptides as positive dataset obtained from DRAMP 
## pos_test_ds1 = 291 antifungal peptides as positive dataset obtained from DRAMP 
## main_antifp_pos_train_ds2 = Dataset contains 1168 antifungal peptides as positive dataset obtained from DRAMP
## amp_otherthan_antifungal _test_ds1 =  291 antimicrobial peptides other than antifungal 
## amp_otherthan_antifungal _train_ds1 = 1168 antimicrobial peptides other than antifungal as negative dataset obtained from DRAMP database.
print("--- Antifp pos_train_ds3  ----")

filenames = ["pos_train_ds1","pos_test_ds1", "amp_otherthan_antifungal _test_ds1",
            "amp_otherthan_antifungal _train_ds1","main_antifp_pos_train_ds2"]
df = pd.DataFrame()

# loop 
for filename in filenames:
    input_path = os.path.join(root_path + "/antifp/", filename+".txt")
    _df = pd.read_csv(input_path, header=None, names=['Sequence'], sep="\t")
    # create new ID
    _df['ID'] = _df.index 
    _df['ID'] = _df['ID'].apply(lambda row: str(row)+"_"+filename)
    df = df.append(_df, ignore_index=True)

# Group duplication
grouped_df = df.groupby('Sequence')
grouped_lists = grouped_df['ID'].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()
filename = "Antifp"
out_path = os.path.join(root_output_path, filename+".fasta")
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

filename = "adam"
input_path = os.path.join(root_path, filename+".txt")
out_path = os.path.join(root_output_path, filename+".fasta")
df = pd.read_csv(input_path, header=0, sep="\t")

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
print("--- End of ADAM ----")

## df[~df["ID"].str.contains("ADAM")]
## df[~df["Sequence"].str.contains("(")]
# %%
# AntiTbPred
print("--- antitbpred  ----")

filename = "antitbpred_pos"
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
        if " " in seq:
           print(seq)
           seq = seq.replace(" ", "")
        file.write(seq + '\n' )
print("--- End of antitbpred ----")

# %%
## YADAMP
print("--- YADAMP ----")
filename = "YADAMP.data"
input_path = os.path.join(root_path, filename+".txt")
out_path = os.path.join(root_output_path, filename+".fasta")

df = pd.read_json(input_path)
# split column in to multiple columns by comma 
df = pd.concat([df['aaData'], df['aaData'].apply(pd.Series)],axis=1)

df = df.drop(['aaData'], axis=1)

# Group duplication  , select on column 6 (Bioactive sequence)
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
        if " " in seq:
           print(seq)
           seq = seq.replace(" ", "")

        # print(seq)
        file.write(seq + '\n')



# %%
## AntiBP2
print("--- antibp2_data ----")
filename = "antibp2"
input_path = os.path.join(root_path + "/antibp2/", filename+".tsv")
out_path = os.path.join(root_output_path, filename+".fasta")
# Load File
df = pd.read_csv(input_path, header=[0], delimiter='\t')

# Group duplication
df['AVP_ID'] = df['AVP_ID'].astype(str)
grouped_df = df.groupby("Sequence")
grouped_lists = grouped_df["AVP_ID"].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        header = createID(row['AVP_ID'], filename)
        # print(header)
        file.write(header + '\n')
        seq = row['Sequence']
        if " " in seq:
            print("Found whitespace:"+seq)
            seq = seq.replace(" ", "")
        # print(seq)
        file.write(seq + '\n')
print("--- End of antibp2_data ----")




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
# Remove non coding 
amp_DB_fasta = os.path.join(DIR,'one_fasta_file.fasta')
PROTIEN_CODE = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y"
allowed_chars  = set("ACDEFGHIKLMNPQRSTVWY")
with open(amp_DB_fasta) as fasta_file:  # Will close handle cleanly
    identifiers = []
    lengths = []
    sequence = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        identifiers.append(seq_record.id)
        sequence.append(str(seq_record.seq))
        lengths.append(len(seq_record.seq))

_AMPs_df = pd.DataFrame(list(zip(identifiers, sequence, lengths)), columns =['ID', 'Sequence', 'length'])
_AMPs_df['Sequence'] = _AMPs_df['Sequence'].str.upper()
#_AMPs_df = _AMPs_df[_AMPs_df["Sequence"].str.contains(PROTIEN_CODE, regex=True)]
# allow some character onyl
_AMPs_df = _AMPs_df[_AMPs_df['Sequence'].apply(lambda x: set(x).issubset(allowed_chars))]
_AMPs_df = _AMPs_df[_AMPs_df["length"]>=5]
# remove unannotated ID
_AMPs_df = _AMPs_df[_AMPs_df["ID"] != "_"]

# get less than or equal 30 
_AMPs_df_30 = _AMPs_df[_AMPs_df["length"] <= 30]
# get greater than 30
_AMPs_df_gt_30 = _AMPs_df[_AMPs_df["length"] > 30]


#### Save 
with open(os.path.join(DIR,'AMP.lte30.fasta'), 'w') as file:
    for index, row in _AMPs_df_30.iterrows():
        file.write('>'+ row['ID'] + '\n')
        file.write(row['Sequence'] + '\n')
file.close()

with open(os.path.join(DIR,'AMP.gt30.fasta'), 'w') as file:
    for index, row in _AMPs_df_gt_30.iterrows():
        file.write('>'+ row['ID'] + '\n')
        file.write(row['Sequence'] + '\n')
file.close()

with open(os.path.join(DIR,'AMP.26DBs.fasta'), 'w') as file:
    for index, row in _AMPs_df.iterrows():
        file.write('>'+ row['ID'] + '\n')
        file.write(row['Sequence'] + '\n')
file.close()





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
