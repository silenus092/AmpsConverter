# AMPsConverter
Tools for fasta conversion that supports the following AMP Databases
1. DBAASP Database of Antimicrobial Activity and Structure of Peptides
2. AVPdb | database of experimentally validated antiviral peptides
3. HIPdb | A Database of Experimentally Validated HIV Inhibiting Peptide
4. Hemolytik Database of Hemolytic activity of AMPs.
5. BaAMPs | Biofilm-active AMPs database: antimicrobial peptides (AMPs) specifically tested against microbial biofilms.
6. MilkAMP A database for milk AMPs (a comprehensive database of antimicrobial peptides of dairy origin)
7. DADP: Database of Anuran Defense Peptides
8. CancerPPD Database of Anticancer peptides and proteins
9. AntiCP Consists of 225 antimicrobial peptides with anticancer properties.
10. LAMP2: An update to LAMP database linking antimicrobial peptide.
11. Antimicrobial Peptide Database (APD3, http://aps.unmc.edu/AP/main.php)
12. EnzyBase  Database of enzybiotics (lysins, autolysins, lysozymes, and large bacteriocins) (Not sure)
13. EROP-Moscow The EROP-Moscow oligopeptide database
14. BACTIBASE : Database Dedicated to Bacteriocin (Filter by antimicrobial keyword  )
15. DRAMP 2.0, an updated data repository of antimicrobial peptides

# Collecting and Fitlering
Strategies:
1. Convert lower case to upper case
2. Remove white space between seqeunce or beginning and ending (e.g. "VEWNEMT WMEWEREI ENYTKLIYKILEESQEQ" to "VEWNEMTWMEWEREIENYTKLIYKILEESQEQ")
3. Remove empty lines
4. Remove none protien seqeunce charater or unrelated words  (e.g. "+FRPKFGK-" or "structure given" or "NOT AVAILABLE")
5. Remove unsure sequence. (e.g. "fCYwO-CyLeu-Pen-TKKr" or "βAlaFTicO")

## **DBAASP**

### <u>How To download original source</u>
1. Visit https://dbaasp.org/search
2. Click "Download CSV" at the bottom of page

### <u>How To convert originial source to FASTA format</u>
```python
root_path = "/mnt/c/works/RKI/AMP_DB/AMPs"
root_output_path = "/mnt/c/works/RKI/AMP_DB/fasta"

def createID(ID, source, seq_type=None, description=None):
    header = str(">" + ID + "|" + source)
    if seq_type is not None:
        header = header + "|" + seq_type
    if description is not None:
        header = header + "|" + description
    return header 

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
```
Note:
 1. lower case?


## **AVPdb**

### <u>How To download original source</u>
1. Visit http://crdd.osdd.net/servers/avpdb/download_data.php

### <u>How To convert originial source to FASTA format</u>
```python
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
```

## **HIPdb**

### <u>How To download original source</u>
1. Visit http://crdd.osdd.net/servers/hipdb/download_data.php

### <u>How To convert originial source to FASTA format</u>
```python
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
```
Note:
 1. file input format (utf-8, utf-16)


## **Hemolytik**	

### <u>How To download original source</u>
1. Download from  http://crdd.osdd.net/raghava/hemolytik/allsequences.txt

### <u>How To convert originial source to FASTA format</u>
```python
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
```
Note:
 1. missing char? `X-NH-GVATQLLAAYILLFDEYNEKKASAQKDILIKVL`
 2. lower case?
 3. βAlaFTicOicβAlaKTicOicβAlaFTicOicβAlaKTicKKKK ?

## **BaAmps**	
### <u>How To download original source</u>
1. Access to this website http://www.baamps.it/browse and select ALL in dropdown, then click export.

### <u>How To convert originial source to FASTA format</u>
```python
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
```
Note:
 1. lower case?

## **MilkAMP**
### <u>How To</u>
1. Download from http://milkampdb.org/entrieslist.php?export=csv

### <u>How To convert originial source to FASTA format</u>
```python
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
        header = createID(row['id'], filename)
        # print(header)
        file.write(header + '\n')
        seq = row['Sequence']
        # print(seq)
        file.write(seq + '\n')

print("--- End of milkampdb ----")
```
Note:
 1. different format? `APRKNVRWCTISQPEW/CIRA;`

## **DADP**
DADP: Database of Anuran Defense Peptides

### <u>How To</u>
Access to the below link and save it as a json file
```
http://split4.pmfst.hr/dadp/?a=listAjax&classId=&search=undefined&limit=undefined
```
### <u>How To convert originial source to FASTA format</u>
```python
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
```


## **CancerPPD**
 Database of Anticancer peptides and proteins

### <u>How To download original source</u>
```bash
curl http://crdd.osdd.net/raghava/cancerppd/natural_seq_dwn/l_natural.txt -o CancerPPD_NaturalAminoAcids_l_natural.txt

curl http://crdd.osdd.net/raghava/cancerppd/natural_seq_dwn/d_natural.txt -o CancerPPD_NaturalAminoAcids_d_natural.txt

curl http://crdd.osdd.net/raghava/cancerppd/non-natural_seq_dwn/l_non-natural.txt -o CancerPPD_l_non-natural.txt

curl http://crdd.osdd.net/raghava/cancerppd/non-natural_seq_dwn/d_non-natural.txt -o CancerPPD_d_non-natural.txt

curl http://crdd.osdd.net/raghava/cancerppd/natural_seq_dwn/mix_natural.txt  -o CancerPPD_mix_natural.txt

curl http://crdd.osdd.net/raghava/cancerppd/non-natural_seq_dwn/mix_non-natural.txt  -o CancerPPD_mix_non-natural.txt
```
### <u>How To convert originial source to FASTA format</u>
```python
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
```


## **AntiCP**
Consists of 225 antimicrobial peptides with anticancer properties.

### <u>How To download original source</u>
1. Access to the below URL and then looks for "Main Dataset"
```
https://webs.iiitd.edu.in/raghava/anticp/datasets.php
```
2. Then click "P" under Donwload column which will redirect to  https://webs.iiitd.edu.in/raghava/anticp/datasets/225_amp_pos.txt

Note:
1. They also provide alternate datasets (e.g., potential non-anticancer peptides , without
anticancer activities)

### <u>How To convert originial source to FASTA format</u>
```python
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
```
## **Antimicrobial Peptide Database (APD3)**

### <u>How To download original source</u>
```
# http://aps.unmc.edu/AP/downloads.php
wget http://aps.unmc.edu/AP/APD3_update2020_release.fasta
```

### <u>How To convert originial source to FASTA format</u>
1. It already offers FASTA format

## **EROP-Moscow**

### <u>How To download original source</u>

1. Go to thise page http://erop.inbi.ras.ru/query-erop.php?
2. Under functional class, select antimicrobial 

or Download from the link below  
```
http://erop.inbi.ras.ru/result1.php?EROP_NMB_K=&PEP_NAME_K=&FAM_NAME_K=&ALL_KAR__K=&ALL_KGD__K=&ALL_PHYL_K=&ALL_B_CL_K=&Organism=&SPECIES__K=&ALL_TISS_K=&SEQ_1____K=&AAR_SUM__K1=&AAR_SUM__K2=&M_W______K1=&M_W______K2=&PI_______K1=&PI_______K2=&POSIT_CH_K1=&POSIT_CH_K2=&NEGAT_CH_K1=&NEGAT_CH_K2=&NET_CH___K1=&NET_CH___K2=&FUNC_CL__K=antimicrobial&FUNCTION_K=&SEQ_REFA_V=&SEQ_REFT_V=&SEQ_REFJ_V=&YEAR_SEQ_V1=&YEAR_SEQ_V2=&COUNTRY__V=&page_mode=Download_Fasta 
```

### <u>How To convert originial source to FASTA format</u>
1. It already offers FASTA format

## **LAMP2**

### <u>How To download original source</u>
1. Download from http://biotechlab.fudan.edu.cn/database/lamp/db/lamp2.fasta

### <u>How To convert originial source to FASTA format</u>
1. It already offers FASTA format

## **BACTIBASE**

### <u>How To download original source</u>
1. Access to  http://bactibase.hammamilab.org/bacteriocinslist.php, and then input the condition on the search bar (e.g. antimicrobial)
2. Select show "all" at the bottom of the page.
3. Click "FASTA" icon to download 

### <u>How To convert originial source to FASTA format</u>
1. It already offers FASTA format

## **EnzyBase**

### <u>How To download original source</u>
1. Downlaod from http://biotechlab.fudan.edu.cn/database/enzybase/db/enzy2.fasta

### <u>How To convert originial source to FASTA format</u>
1. It already offers FASTA format

## **DRAMP (2.0)**

### <u>How To</u>
1. Vist the http://dramp.cpu-bioinfor.org/downloads/
2. Click download link

### <u>How To convert originial source to FASTA format</u>
```python
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

print("--- DRAMP ----")

for file_name in file_names:
    print(file_name)
    converterForDRAMP(file_name)

print("--- End of DRAMP ----")
```
