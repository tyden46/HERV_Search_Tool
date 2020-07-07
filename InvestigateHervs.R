library(Rgb)
library(biomaRt)
library(rentrez)

#You can download HERV annotation at:
#https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg38.v2/transcripts.gtf 
HERVtranscripts=read.gtf("transcripts.gtf") #Matthew's HERV annotation
res.sig=read.csv("AlzheimersHervs-PRJEB28518.csv", row.names = 1)
#Step 1 : Find Genomic Coordinates for HERVs
tableExists=FALSE #Create table first, then rbind to it
for(x in rownames(res.sig)){ #Iterate through sig HERVs
  idx=which(HERVtranscripts$gene_id==x)
  if(length(idx)>1){ #If there are multiple loci for a herv, pick the first
    thisIndex=idx[[1]][1]
  }else{
    thisIndex=idx
  }
  #9th col = gene_id
  #1st col = chromosome
  #4th col = start
  #5th col = end
  #By modifying the values in cols 4&5 we could look upstream
  #Or downstream of a HERV location to search for genes
  if(tableExists==FALSE){
    myTable=as.data.frame(HERVtranscripts[thisIndex, c(9,1,4,5)])
    tableExists=TRUE
  }else{
    myTable=as.data.frame(rbind(myTable,HERVtranscripts[thisIndex, c(9,1,4,5)]))
  }
}

#Search for genes in HERV locations using biomaRt
myMart <- useMart(biomart = "ensembl", 
                  host="www.ensembl.org",
                  dataset = "hsapiens_gene_ensembl") 
attributes <- c("ensembl_gene_id","hgnc_symbol") #We could fetch more attributes
filters <- c("chromosome_name","start","end")
listOfGenes=c()
listOfValues=c()
#Here we construct our query list
for(x in 1:length(myTable$gene_id)){
  thisChromosome=myTable[x,2] #Will fetch "chrX" "chr17" etc
  #We need to remove "chr" substring
  chromosomeNumber=substr(thisChromosome, 4, nchar(thisChromosome))
  startNum=myTable[x,3]
  endNum=myTable[x,4]
  #We create a vector of search parameters
  values <- (c(chromosome=chromosomeNumber,
               start=startNum,
               end=endNum))
  #We add to our list of genomic queries
  listOfValues=c(listOfValues, entry=list(values))
}
#Now that we have a query list we can search using biomaRt
tableExists=FALSE
#for(x in 1:length(row.names(myTable))){
for(x in 1:10){
  print(paste("Searching for gene number", x))#Acts as a counter so you can see progress, this loop takes a while
  #We pass in each query to get the attributes we want
  foundGenes=FALSE
  extendFactor=0
  while(foundGenes==FALSE && extendFactor<10000){
    chrom=as.character(listOfValues[[x]][1])
    startPos=start = as.character(as.numeric(as.character(listOfValues[[x]][2]))-extendFactor)
    endPos=start = as.character(as.numeric(as.character(listOfValues[[x]][3]))+extendFactor)
    theseGenes <- getBM(attributes=attributes,
                        filters=filters,
                        values=list(chrom,
                                    startPos,
                                    endPos),
                        mart=myMart)
    if(length(theseGenes$ensembl_gene_id)>0){
      foundGenes=TRUE
    }else{
      print("Extending search space")
      extendFactor=extendFactor+500
    }
  }
  if(tableExists){#If table exists, add entry
    geneTable=rbind(geneTable, theseGenes)    
  }else{#If table doesn't exist yet, make it
    geneTable=theseGenes
    tableExists=TRUE
  }
  foundGenes=FALSE
}
#Now that we have a list of genes we can search pubmed
for(x in geneTable$hgnc_symbol){
  if(x!= "" && !is.na(x)){#If the gene has a symbol
    #We could cast a wider net here by searching the abstracts as well as title
    #We could search for "brain" or "neur*"
    myTerm=paste("(", x, "[Title/Abstract]) AND ((Alzheimer's[Title/Abstract]) OR (neur*[Title/Abstract]) OR (brain[Title/Abstract]))", sep="")
    #myTerm=paste("(", x, "[Title]) AND (Alzheimer's[Title])", sep="")
    search <- entrez_search(db="pubmed",
                            term=myTerm)
    #If results are found, print how many papers match the search
    if(length(search$ids)>0){
      print(paste("Found", length(search$ids),
                  "results for", x, sep=" "))
    }
  }
}