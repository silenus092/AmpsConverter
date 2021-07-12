# TO generate likely AMP set

## MycoCosm

```bash
#!/bin/bash
### AMPs  > 30

# Search against
mmseqs easy-search /mnt/vdb/data/26dbs/AMP.gt30.fasta /mnt/vdb/DBs/fungi/DB/targetDB /mnt/vdb/DBs/fungi/amp-fungi.gt30.nr /mnt/vdb/DBs/fungi/DB/tmp  -s 7.5 -c 1  --split-memory-limit  20G
# mmseqs easy-search /mnt/vdb/data/26dbs/AMP.gt30.fasta /mnt/vdb/DBs/fungi/DB/targetDB /mnt/vdb/DBs/fungi/amp-fungi.gt30.nr /mnt/vdb/DBs/fungi/DB/tmp  -s 7.5 -c 0.9  --split-memory-limit  20G
# mmseqs easy-search /mnt/vdb/data/26dbs/AMP.gt30.fasta /mnt/vdb/DBs/fungi/DB/targetDB /mnt/vdb/DBs/fungi/amp-fungi.gt30.nr /mnt/vdb/DBs/fungi/DB/tmp  -s 7.5 -c 0.95  --split-memory-limit  20G


# fitler by at least 95% similarity and then save into 
awk '{ if ($3 >= 0.985) { print $2 } }' /mnt/vdb/DBs/fungi/amp-fungi.gt30.nr > /mnt/vdb/DBs/fungi/fungi.gt30.sim95.txt

# get uniqueID
sort /mnt/vdb/DBs/fungi/fungi.gt30.sim95.txt | uniq > /mnt/vdb/DBs/fungi/fungi.gt30.sim95.uniq.txt

# extract from list of ID  
seqkit  grep -j 24 -f /mnt/vdb/DBs/fungi/fungi.gt30.sim95.uniq.txt /mnt/vdb/DBs/2020-04-15_all_fungi_proteins.faa.gz  -o /mnt/vdb/DBs/fungi/amp_like.gt30.sim95.fungi.fa
```

```bash
#!/bin/bash
### AMPs  <= 30
## Fungi
mmseqs  easy-search /mnt/vdb/data/26dbs/AMP.lte30.fasta /mnt/vdb/DBs/fungi/DB/targetDB /mnt/vdb/DBs/fungi/amp-fungi.lte30.nr /mnt/vdb/DBs/fungi/DB/tmp  -s 7.5 -c 0.95 --split-memory-limit 20G
#mmseqs  easy-search /mnt/vdb/data/26dbs/AMP.lte30.fasta /mnt/vdb/DBs/fungi/DB/targetDB /mnt/vdb/DBs/fungi/amp-fungi.lte30.nr /mnt/vdb/DBs/fungi/DB/tmp  -s 7.5 -c 0.1 --split-memory-limit 20G

# fitler by at least 80% similarity and then save into 
awk '{ if ($3 >= 0.98) { print $2 } }' /mnt/vdb/DBs/fungi/amp-fungi.lte30.nr > /mnt/vdb/DBs/fungi/fungi.lte30.sim98.txt

# get uniqueID
sort /mnt/vdb/DBs/fungi/fungi.lte30.sim98.txt | uniq > /mnt/vdb/DBs/fungi/fungi.lte30.sim98.uniq.txt

# extract from list of ID  
seqkit  grep -j 24 -f /mnt/vdb/DBs/fungi/fungi.lte30.sim98.uniq.txt /mnt/vdb/DBs/2020-04-15_all_fungi_proteins.faa.gz  -o /mnt/vdb/DBs/fungi/amp_like.lte30.sim98.fungi.fa
```

## MetaClust
running at EBI big mem node YODA *
```bash
module load singularity/3.5.0
# screen on mhoelzer@hh-yoda-08-01

#bsub -P bigmem -n 24 -M 1500 -R rusage[mem=1500] -Is bash
# hh-yoda-bm-01
bsub -P bigmem -n 24 -M 900000 -R rusage[mem=900000] -Is bash

cd /hps/nobackup2/metagenomics/mhoelzer/data/2021-02-20_amp_metaclust_extension

wget https://metaclust.mmseqs.com/current_release/metaclust_nr.fasta.gz
gzip -d metaclust_nr.fasta.gz

# 22DBs search against source  (conda install -c conda-forge -c bioconda mmseqs2) if required
singularity pull --name mmseqs.img docker://nanozoo/mmseqs2:11.e1a1c--55acb62

for TYPE in gt30 lte30; do for C in 1 0.95 0.9; do

# Seach 
singularity run mmseqs.img mmseqs easy-search -s 7.5 -c $C AMP.$TYPE.fasta metaclust_nr.fasta amp-metaclust.$TYPE.$C.nr tmp

# Fitler at least 98.5% similarity 
awk '{ if ($3 >= 0.98) { print $2 } }' amp-metaclust.$TYPE.$C.nr > metaclust.$TYPE.$C.txt

sort metaclust.$TYPE.$C.txt | uniq > metaclust.$TYPE.$C.uniq.txt

# extract from list of ID  
seqkit grep -j 24 -f metaclust.$TYPE.$C.uniq.txt metaclust_nr.fasta -o metaclust.$TYPE.$C.fa

done; done
```

## MGnify

```bash
#!/bin/bash
### AMPs  > 30
# Search against
mmseqs  easy-search /mnt/vdb/data/26dbs/AMP.gt30.fasta /mnt/vdb/DBs/mgnify/DB/targetDB /mnt/vdb/DBs/mgnify/amp-mgnify.gt30.nr /mnt/vdb/DBs/mgnify/DB/tmp  -s 7.5 -c 1 --split-memory-limit  30G
#mmseqs  easy-search /mnt/vdb/data/26dbs/AMP.gt30.fasta /mnt/vdb/DBs/mgnify/DB/targetDB /mnt/vdb/DBs/mgnify/amp-mgnify.gt30.nr /mnt/vdb/DBs/mgnify/DB/tmp  -s 7.5 -c 0.95 --split-memory-limit  30G
#mmseqs  easy-search /mnt/vdb/data/26dbs/AMP.gt30.fasta /mnt/vdb/DBs/mgnify/DB/targetDB /mnt/vdb/DBs/mgnify/amp-mgnify.gt30.nr /mnt/vdb/DBs/mgnify/DB/tmp  -s 7.5 -c 0.9 --split-memory-limit  30G

# fitler by at least 80% similarity and then save into 
awk '{ if ($3 >= 0.985) { print $2 } }' /mnt/vdb/DBs/mgnify/amp-mgnify.gt30.nr > /mnt/vdb/DBs/mgnify/mgnify.gt30.sim95.txt

#awk '{ if ($3 >= 0.98) { print $2 } }' /mnt/vdb/DBs/mgnify/c95/amp-mgnify.gt30.nr > /mnt/vdb/DBs/mgnify/c95/mgnify.gt30.txt
#awk '{ if ($3 >= 0.98) { print $2 } }' /mnt/vdb/DBs/mgnify/c90/amp-mgnify.gt30.nr > /mnt/vdb/DBs/mgnify/c90/mgnify.gt30.txt


# get uniqueID
sort /mnt/vdb/DBs/mgnify/mgnify.gt30.sim95.txt | uniq > /mnt/vdb/DBs/mgnify/mgnify.gt30.sim95.uniq.txt
#sort /mnt/vdb/DBs/mgnify/c95/mgnify.gt30.txt | uniq > /mnt/vdb/DBs/mgnify/c95/mgnify.gt30.uniq.txt
#sort /mnt/vdb/DBs/mgnify/c1/mgnify.gt30.txt | uniq > /mnt/vdb/DBs/mgnify/c1/mgnify.gt30.uniq.txt

# extract from list of ID  
seqkit  grep -j 24 -f /mnt/vdb/DBs/mgnify/mgnify.gt30.sim95.uniq.txt /mnt/vdb/DBs/mgy_clusters_20190531.fa.gz -o /mnt/vdb/DBs/mgnify/amp_like.gt30.sim95.mgnify.fa
#seqkit  grep -j 24 -f /mnt/vdb/DBs/mgnify/c95/mgnify.gt30.uniq.txt /mnt/vdb/DBs/mgy_clusters_20190531.fa.gz -o /mnt/vdb/DBs/mgnify/c95/amp_like.gt30.mgnify.fa
#seqkit  grep -j 24 -f /mnt/vdb/DBs/mgnify/c90/mgnify.gt30.uniq.txt /mnt/vdb/DBs/mgy_clusters_20190531.fa.gz -o /mnt/vdb/DBs/mgnify/c90/amp_like.gt30.mgnify.fa
```


```bash
#!/bin/bash
### AMPs  <= 30
mmseqs  easy-search  /mnt/vdb/data/26dbs/AMP.lte30.fasta /mnt/vdb/DBs/mgnify/DB/targetDB /mnt/vdb/DBs/mgnify/amp-mgnify.lte30.nr /mnt/vdb/DBs/mgnify/DB/tmp -s 7.5 -c 0.95 --split-memory-limit 30G

# fitler by at least 80% similarity and then save into 

awk '{ if ($3 >= 0.98) { print $2 } }' /mnt/vdb/DBs/mgnify/amp-mgnify.lte30.nr > /mnt/vdb/DBs/mgnify/mgnify.lte30.sim98.txt
# awk '{ if ($3 >= 0.90) { print $2 } }' /mnt/vdb/DBs/mgnify/c95/amp-mgnify.lte30.nr > /mnt/vdb/DBs/mgnify/c95/mgnify.lte30.txt


# get uniqueID

sort /mnt/vdb/DBs/mgnify/mgnify.lte30.sim98.txt | uniq > /mnt/vdb/DBs/mgnify/mgnify.lte30.sim98.uniq.txt
# sort /mnt/vdb/DBs/mgnify/c95/mgnify.lte30.txt | uniq > /mnt/vdb/DBs/mgnify/c95/mgnify.lte30.uniq.txt


# extract from list of ID  

seqkit  grep -j 24 -f /mnt/vdb/DBs/mgnify/mgnify.lte30.sim98.uniq.txt /mnt/vdb/DBs/mgy_clusters_20190531.fa.gz -o /mnt/vdb/DBs/mgnify/amp_like.lte30.sim98.mgnify.fa
# seqkit  grep -j 24 -f /mnt/vdb/DBs/mgnify/c95/mgnify.lte30.uniq.txt /mnt/vdb/DBs/mgy_clusters_20190531.fa.gz -o /mnt/vdb/DBs/mgnify/c95/amp_like.lte30.mgnify.fa
```
# Merging

* Merge all external sources and 26.AMP.fasta together 
* Clean and remove ambiguous amino acid coding
* deduplicate records

# Cluster

Set A
```
mmseqs createdb amp.ready.fasta DB

mmseqs cluster -c 1 --min-seq-id 0.4 -s 7.5 DB DB_clu tmp

mmseqs createtsv DB DB DB_clu DB_clu.tsv

mmseqs createseqfiledb DB DB_clu DB_clu_seq

mmseqs result2flat DB DB DB_clu_seq final_AMP_set_A.fasta
```

Set B
```
mmseqs createdb amp.ready.fasta DB

mmseqs cluster -c 0.4 --min-seq-id 0.4 -s 7.5 DB DB_clu tmp

mmseqs createtsv DB DB DB_clu DB_clu.tsv

mmseqs createseqfiledb DB DB_clu DB_clu_seq

mmseqs result2flat DB DB DB_clu_seq final_AMP_set_B.fasta
```