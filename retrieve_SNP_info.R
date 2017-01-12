#first, load this function in R:
#R
#source("retrieve_SNP_info.R")

#then, to run this function:
#retireve_SNP_info(snp_id)
#where snp_id is the SNP name or rsid, like:
#retrieve_SNP_info("rs10828658")

#if you want to store the results in an object
#my_result <- retrieve_SNP_info("rs10828658")

retrieve_SNP_info <- function (snp_id) {

	#load biomaRt library
	library(biomaRt)

	print("Loading SNP database...")

	#select database
	snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
	ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

	#selects what attributes to show, removing some buggy attributes that cause errors
	snp_attributes = listAttributes(snp_mart)$name[-c(65,66,67)]

	#results will be stored in this object
	snp_table <- 0

	#here are the indexes of some key info from snp_attributes
	#37 -> ensembl_gene_stable_id
	#66 -> chrom_start
	#67 -> chr_name
	#11 -> allele_1
	#12 -> minor_allele

	print("Loading SNP info...")
	
	#if you want to see all SNP info, uncomment this line bellow, and comment the next for loop
	#for(i in 1:length(snp_attributes)){

	#or if you want to just check just some key info, leave the line below uncomented, and leave the line above commented
	for(i in c(37,66,67,11,12)){

		#print(i)
		#print(snp_attributes[i])
		tmp <- getBM(attributes=snp_attributes[i],filters="snp_filter", values=snp_id, mart=snp_mart)


		#test if result is empty
		if(nrow(tmp)!=0){

			if(snp_table!=0){
				#copy info to snp_table, if information is not empty
				snp_table <- cbind(snp_table, tmp)				
			} else {
				snp_table <- tmp
			}

			#print(tmp)
		}	

	}

	print("Loading Ensembl gene database...")

	#load ensambl gene database for homo sapiens
	ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

	#create a variable with a list of gene attributes
	gene_attributes <- listAttributes(ensembl)$name
	
	#here are the indexes of some key info from gene_attributes
	#1 -> ensembl_gene_id
	#1382 -> Chromosome/scaffold name
	#1383 -> Gene Start (bp)
	#1384 -> Gene End (bp)

	gene_table <- 0

	print("Loading gene info...")

	#if you want to see all gene info, uncomment this line bellow, and comment the next for loop
	#for(i in 1:length(gene_attributes)){

	#get the chromosome name, start, and end position for this gene
	#or if you want to just check just some key info, leave the line below uncomented, and leave the line above commented
	for(i in c(1,1382,1383,1384)){

		#print(i)
		#print(gene_attributes[i])
		tmp <- getBM(attributes = gene_attributes[i], values = snp_table[1,"ensembl_gene_stable_id"], filters=c(gene_attributes[1]), mart = ensembl)


		#test if result is empty
		if(nrow(tmp)!=0){

			if(gene_table!=0){
				#copy info to gene_table, if information is not empty
				gene_table <- cbind(gene_table, tmp)				
			} else {
				gene_table <- tmp
			}

			#print(tmp)
		}	

	}
	
	
	print("Processing gene sequence variations...")

	#retrieve the gene sequence 
	SNP_gene_sequence = getSequence(id=snp_table[1,"ensembl_gene_stable_id"], type=gene_attributes[1], seqType="gene_exon_intron", mart = ensembl)

	#calculate the SNP location inside the gene
	snp_gene_pos <- snp_table[1, "chrom_start"]-gene_table[1,"start_position"]+1


	#if you want to check whether the SNP is located in the correct position
	if(substr(SNP_gene_sequence[1], snp_gene_pos, snp_gene_pos) != snp_table[1, "allele_1"]){
		print("Error gene location doesn't match to allele_1!")
		return(1)
	}

	#generate 2 copies of the gene sequence, one will keep the ancestral allele, the other will have the minor allele
	SNP_gene_sequence_var1 <- as.character(SNP_gene_sequence[1])
	SNP_gene_sequence_var2 <- as.character(SNP_gene_sequence[1])

	#replace the base in the SNP location by "minor_allele" base
	substr(SNP_gene_sequence_var2, snp_gene_pos, snp_gene_pos) <- snp_table[1,"minor_allele"]


	#store results into a data frame
	result1 <- data.frame(SNP=snp_id, VARIATION="var1", SEQUENCE=SNP_gene_sequence_var1)
	result2 <- data.frame(SNP=snp_id, VARIATION="var2", SEQUENCE=SNP_gene_sequence_var2)

	#merge results
	result <- rbind(result1,result2)

	print("Done!")

	return(result)
}
