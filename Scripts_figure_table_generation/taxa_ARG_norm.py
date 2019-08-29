#!/bin/biopython 
import pandas as pd 
import sys 
#import ARG_kraken
data=pd.read_csv(sys.argv[1],header=None,sep="\t")
dict_reads_taxa=dict(zip(data[0].tolist(),data[1].tolist()))
#import ARG table 
ARG_table=pd.read_csv(sys.argv[2],header=None,sep="\t")
reads_taxa=data[0].tolist()
ARG_table_filter=ARG_table.loc[ARG_table[1].isin(reads_taxa)]
reads_gene_dict=dict(zip(ARG_table_filter[1].tolist(),ARG_table_filter[2].tolist()))
from Bio import SeqIO
fasta_sequences = SeqIO.parse(open('megares_database_v1.01.fasta'),'fasta')
gene=[]
seq_len=[]
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    gene.append(name)
    seq_nt_number=len(sequence)
    seq_len.append(seq_nt_number)
    
gene_name=[item.split("|")[-1] for item in gene]
ARG_gene_dict=dict(zip(gene_name,seq_len))

reads=[]
length=[]
for key, value in reads_gene_dict.items() :
    reads.append(key)
    length.append(ARG_gene_dict[value])
reads_length_dict=dict(zip(reads,length))
taxa=[]
length_to_calculate=[]
for key, value in reads_length_dict.items():
    taxa.append(dict_reads_taxa[key])
    length_to_calculate.append(value)
d=list(zip(taxa, length_to_calculate))
df= pd.DataFrame(d, columns =['taxa', 'length']) 
#output taxa and gene length 
df.to_csv(sys.argv[3],index=None,sep="\t")


