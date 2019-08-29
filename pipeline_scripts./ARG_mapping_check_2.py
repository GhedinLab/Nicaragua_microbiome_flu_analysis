import pandas as pd 
import sys
#read in the alignment data 
data=pd.read_csv(sys.argv[1],sep="\t",header=None,names = ["Ref", "Start", "End", "Reads","Score","Strand"])
#set apart the data needs SNP confirmation 
data_1=data[data["Ref"].str.contains('RequiresSNPConfirmation')]
data_2=data[data["Ref"].str.contains('RequiresSNPConfirmation') == False]
#only looking at the data with confirmed antibiotic resistance genes
Ref=data_2["Ref"].tolist()
Ref_gene=[x.split("|")[-1] for x in Ref]
Ref_gene_class=[x.split("|")[-3] for x in Ref]
#get the table for antibiotic resistance genes and the classes of antibiotics the genes confer resistance to
Ref_table=pd.DataFrame({'Gene':Ref_gene,'class':Ref_gene_class})
#remove the duplicate entities
Ref_data=Ref_table.drop_duplicates(subset=None, keep='first', inplace=False)
#create a list of reference genes that will be used to get the alignment information later 
Ref_gene_list=list(set(Ref_gene))
#get the coordinates of the alignment, ref gene and reads aligned to the reference 
data_Resistance_gene_1=pd.DataFrame({'Ref':Ref_gene,'Start':data_2["Start"],'End':data_2['End'],'Reads':data_2["Reads"]})
#remove the duplicate alignments
data_Resistance_gene=data_Resistance_gene_1.drop_duplicates(subset=None, keep='first', inplace=False)
#create a list for the start coordinate, end coordinate, the distance between start and end and the number of reads mapped to each reference gene
#The homologous genes will be count as one gene. 
Start_list=list()
End_list=list()
SE_list=list()
Reads_number=list()
#get the information for the ref genes, start and end coordinates of alignments 
for item in Ref_gene_list:
    #for each reference gene in the table, get the coordinate of the alignments
    Table=data_Resistance_gene[data_Resistance_gene['Ref']== item]
    S=Table["Start"].tolist()
    S=[int(x) for x in S]
    S.sort()
    #for each reference gene, get the first coordinate the reads mapped to 
    Start_list.append(S[0])
    E=Table["End"].tolist()
    E=[int(x) for x in E]
    #for each reference gene, get the last coordinate on the gene that the reads mapped 
    E.sort()
    End_list.append(E[-1])
    #check the mapping distance, is the last coordinate minus the first coordinate. This mapping distance will be used to filter for the ARG genes detected. Only the genes being hitted at least two regions will be keeped. 
    SE_list.append(E[-1]-S[0])
    #count the reads hitting the gene by counting the number of the start coordinate
    Reads_number.append(len(S))
#output the table 
#output=pd.DataFrame({'Ref':Ref_gene_list,'Mapping':SE_list,'Start':Start_list,'End':End_list,'Reads_number':Reads_number})
#output.to_csv(sys.argv[2],sep="\t",index=False)
#output the table for the antibiotic resistance genes and the corresponding classes
#Ref_data.to_csv(sys.argv[3],sep="\t",index=False)
#The gene table
data_Resistance_gene.to_csv(sys.argv[2],sep="\t",index=False)

