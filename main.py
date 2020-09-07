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
"""
# %%
import pandas as pd
import os 

root_path = "/mnt/c/works/RKI/AMP_DB/AMPs"
root_output_path = "/mnt/c/works/RKI/AMP_DB/fasta"


def createID(ID, source, seq_type=None, description=None):
    header = str(">" + ID + "|" + source)
    if seq_type is not None:
        header = header + "|" + seq_type
    if description is not None:
        header = header + "|" + description
    return header 


# %%

# DRAMP
def converterForDRAMP(filename):
    input_path = os.path.join(root_path, filename+".xlsx")
    out_path = os.path.join(root_output_path, filename+".fasta")
    # Load File
    df = pd.read_excel(input_path, header=[0])
    df = groupDRAMPDuplication(df)
    with open(out_path, 'w') as file:
        for index, row in df.iterrows():
            header = createID(row['DRAMP_ID'], filename)
            #print(header)
            file.write(header + '\n')
            seq = row['Sequence']
            #print(seq)
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
def converterForCancerPPD(filename):
    input_path = os.path.join(root_path, filename+".txt")
    out_path = os.path.join(root_output_path, filename+".fasta")
    # Load File
    df = pd.read_csv(input_path, sep='\t', header=[0])
    df = groupCancerPPDDuplication(df)
    with open(out_path, 'w') as file:
        for index, row in df.iterrows():
            header = createID(str(row['id']), filename)
            #print(header)
            file.write(header + '\n')
            seq = row['Sequence'].upper()
            #print(seq)
            file.write(seq + '\n')

def groupCancerPPDDuplication(df):
    df['id'] = df['id'].astype(str)
    grouped_df = df.groupby("Sequence")
    grouped_lists = grouped_df["id"].agg(lambda column: ",".join(column))
    grouped_lists = grouped_lists.reset_index(name="id")
    return grouped_lists

file_names=["CancerPPD_d_natural",
            "CancerPPD_d_non-natural",
            "CancerPPD_l_natural",
            "CancerPPD_l_non-natural",
            "CancerPPD_mix_natural",
            #"CancerPPD_mix_non-natural"
            ]

print("--- CancerPPD ----")

for file_name in file_names:
    print(file_name)
    converterForCancerPPD(file_name)

print("--- End of CancerPPD ----")


# %%
input_path = os.path.join(root_path, "CancerPPD_mix_natural.txt")
df = pd.read_csv(input_path, sep='\t', header=[0])
df['id'] = df['id'].astype(str)
grouped_df = df.groupby("Sequence")
grouped_lists = grouped_df["id"].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index(name="id")



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

print(file_name)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        header = createID(row['Id'], filename)
        # print(header)
        file.write(header + '\n')
        seq = row['Sequence']
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
        seq = row['PeptideSequence']
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

# Group duplication
df['ID'] = df['ID'].astype(str)
grouped_df = df.groupby("SEQ")
grouped_lists = grouped_df["ID"].agg(lambda column: ",".join(column))
grouped_lists = grouped_lists.reset_index()

print(filename)
with open(out_path, 'w') as file:
    for index, row in grouped_lists.iterrows():
        header = createID(row['ID'], filename)
        # print(header)
        file.write(header + '\n')
        seq = row['SEQ']
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
        # print(seq)
        file.write(seq + '\n')
print("--- End of HIPdb_data ----")


# %% 
# dbaasp_peptides
# lower case?
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
        seq = row['SEQUENCE']
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