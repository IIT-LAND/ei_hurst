# genelistOverlap.R
#
# Calculate enrichment odds ratio and p-value from hypergeometric test to
# answer the question of whether genes from one list are enriched in genes
# from another list.
#
# INPUT
#	list1 or list2 = excel file, tab or comma delimited file with gene IDs
#					 assuming each list has a header.
#	backgroundTotal = specify number of genes in background pool
#
# mvlombardo 21.12.2016

genelistOverlap <- function(list1,list2,backgroundTotal, print_result = TRUE, header = FALSE) {
	
	# Read in libraries and set options
	options(stringsAsFactors = FALSE)
	require(readxl)
	require(tools)
  
  if (is.character(list1)){
    # get the file extension of list1
    ext1 = file_ext(list1)
    
    if (is.element(ext1,c("xlsx","xls","txt","csv"))){
      if (ext1=="xlsx" | ext1=="xls") {
        genes1 = read_excel(list1)
      } else if (ext1=="txt") {
        genes1 = read.delim(list1, header = header)
      } else if (ext1=="csv") {
        genes1 = read.csv(list1, header = header)
      }# if (ext1=="xlsx" | ext1=="xls") {
    } else {
      genes1 = data.frame(list1)
    }# if (is.element(ext1,c("xlsx","xls","txt","csv"))){
  } else if (is.data.frame(list1)){
    genes1 = list1
  }# if (is.character(list1)){
  
  if (is.character(list2)){
    # get the file extension of list1
    ext2 = file_ext(list2)
    
    if (is.element(ext2,c("xlsx","xls","txt","csv"))){
      if (ext2=="xlsx" | ext2=="xls") {
        genes2 = read_excel(list2)
      } else if (ext2=="txt") {
        genes2 = read.delim(list2, header = header)
      } else if (ext1=="csv") {
        genes2 = read.csv(list2, header = header)
      }# if (ext1=="xlsx" | ext1=="xls") {
    } else {
      genes2 = data.frame(list2)
    }# if (is.element(ext1,c("xlsx","xls","txt","csv"))){
  } else if (is.data.frame(list2)){
    genes2 = list2
  }# if (is.character(list2)){
  
	# Find overlapping genes
	gene_mask = is.element(genes1[,1],genes2[,1])
	overlapping_genes = genes1[gene_mask,1]
	gene_overlap = sum(gene_mask)
	ngenes1 = length(genes1[,1])
	ngenes2 = length(genes2[,1])

	# Calculate odds ratio
	A = gene_overlap;
	B = ngenes1-gene_overlap
	if (ngenes2==gene_overlap){
	  # add 0.5 to ngenes2 to avoid OR = Inf
	  C = (ngenes2+0.5)-gene_overlap
	} else {
	  C = ngenes2-gene_overlap
	}
	D = backgroundTotal-C
	OR = (A*D)/(B*C)

	# Calculate p-value from hypergeometric test
	hypergeo_p = sum(dhyper(gene_overlap:ngenes2,ngenes1,backgroundTotal-ngenes1,ngenes2))

	# pack into result
	result = vector(mode = "list", length = 1)
	result[[1]]$list1 = list1
	result[[1]]$list2 = list2
	result[[1]]$backgroundTotal = backgroundTotal
	result[[1]]$OR = OR
	result[[1]]$hypergeo_p = hypergeo_p
	result[[1]]$percent_overlap_list1 = gene_overlap/ngenes1
	result[[1]]$percent_overlap_list2 = gene_overlap/ngenes2
	result[[1]]$gene_overlap = gene_overlap
	result[[1]]$ngenes1 = ngenes1
	result[[1]]$ngenes2 = ngenes2
	result[[1]]$overlapping_genes = overlapping_genes

	# print result to the screen and then return result
	if (print_result){
	  print(sprintf("OR = %f, p = %f",OR,hypergeo_p))
	}
	return(result)
} # function genelistOverlap 