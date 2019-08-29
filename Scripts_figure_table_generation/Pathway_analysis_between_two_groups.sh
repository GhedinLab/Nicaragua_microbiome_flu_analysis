#!/bin/bash 

perl /Users/lingdizhang/Desktop/FMAP-master-2/FMAP_comparison.pl gene.table.txt '30082,30112,30191,30291,30292,30295,6019,6205,6471,6628,6667,6777,6929,6957,7047,7098,7331,7637,7965,8033' '30304,30348,30427,30429,30455,30515,3522,6271,6617,7153,7765,30076,8561'  -f 2 -p 0.05 -a 0.05 > gene_partition_ARG_comparions_adjust.txt 

perl /Users/lingdizhang/Desktop/FMAP-master-2/FMAP_pathway.pl gene_partition_ARG_comparions_adjust.txt > gene_partition_pathway_adjust.txt
