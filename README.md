An R script that retrieves the complete gene sequence and its variations using a SNP name or rsid as input. The script is available at Github:

Download the script from:
https://github.com/raqueldias/retrieve_gene_variations_from_SNP/blob/master/retrieve_SNP_info.R

The script requires the biomaRt package, so you gotta install it first in R:

> source("https://bioconductor.org/biocLite.R") 

> biocLite("biomaRt")

After installing the biomaRt package and downloading the script retrieve_SNP_info.R , load and run the R script:

> source("retrieve_SNP_info.R")

> my_result <- retrieve_SNP_info("rs10828658")

Then save the results to a tabular text file:

> write.table(my_result, file="retrieve_gene_result.txt", quote=F, row.names=F, sep="\t")

A full example of result file is available at Github: 

https://raw.githubusercontent.com/raqueldias/retrieve_gene_variations_from_SNP/master/retrieve_result_example.txt

The result contains the complete gene sequence (exons and introns) variations for the gene that belongs to the same region as the SNP rs10828658. 


