

# I downlaoded the data at DNA, RNA and protein levels from CellMiner
# (processed data)
# I removed from the .rda object returned by thsu script 
# the kinase proteomic profile since the naumber of NA values 
# per sample was excessively large.


require(MultiAssayExperiment)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

prep.protein.ms.assay <- function(assayf){
 assayf = assayf
 assayf[assayf==0] = NA
 assay <- assayf[which(apply(assayf[,grep('iBAQ',
  colnames(assayf))],1, function(x) !all(is.na(x)))),]
# assay = as.data.frame(assay)
 assay = as.matrix(assay)
 colnames(assay) <- colnames(assayf)

# feats <- unlist(lapply(lapply(as.character(assay$Fasta.headers),
 feats <- as.character(unlist(lapply(lapply(as.character(assay[,7]),
  function(x) (strsplit(x,"Gene_Symbol=")[[1]][2])), function(x)
  (strsplit(x," ")[[1]][1]))))
 uniquefeats <- names(table(feats))[as.numeric(table(feats))==1]
 assay <- assay[match(uniquefeats,feats),]
 rownames(assay) <- feats[match(uniquefeats,feats)]

 assay <- assay[,grep('iBAQ',colnames(assay))]
 colnames(assay) <- gsub('\\.','\\:',colnames(assay))

 return(as.matrix(assay))
}

# *************************************************** #
# Acquire expression data and cell line-related pdata #
# *************************************************** #

# cell lines meta data
pdataf = read.delim('./nci60_cell_line_metadata.txt',sep='\t',comment.char='#',
 header=T,as.is=T)
# Downloaded from CellMiner on 14/03/2016
rnaf = read.table('./nci60_RNA_Affy_HG_U133_Plus_2.0_RMA/RNA_Affy_HG_U133_Plus_2.0_RMA.txt',sep='\t',header=T,as.is=T)
# Downloaded from CellMiner on 14/03/2016
protein.rpla.f = read.table('./nci60_Protein_Lysate_Array_log2/Protein_Lysate_Array_log2.txt',sep='\t',header=T,as.is=T)
# Downloaded from PMID:23933261
protein.ms.f = as.matrix(read.table('./nci60_Protein_pmid_23933261/pmid.23933261.ST3.nci60.proteomes.txt',sep='\t',header=T,as.is=T))
deep.protein.ms.f = read.table('./nci60_Protein_pmid_23933261/pmid.23933261.ST3.nci60.deep.proteomes.txt',sep='\t',header=T,as.is=T)
kinase.protein.ms.f = read.table('./nci60_Protein_pmid_23933261/pmid.23933261.ST3.nci60.kinomes.txt',sep='\t',header=T,as.is=T)
dnaf = read.table('./nci60_DNA_aCGH_Agilent_44K_gene_summary.txt',
 sep='\t',header=T,as.is=T)

# ************************************************************************** #
# format multi-level and multi-platform expression data and save in R object #
# ************************************************************************** #
 
rna.arraydat = as.matrix(rnaf[,10:69])
rna.arraydat = trim(rna.arraydat)
rownames(rna.arraydat) = rnaf$Gene.name.c
colnames(rna.arraydat) = gsub('\\.','\\:',colnames(rnaf)[10:69])
rna.arraydat = rna.arraydat[match(unique(rownames(rna.arraydat)),
 rownames(rna.arraydat)),]

protein.rpla.arraydat = as.matrix(protein.rpla.f[,10:69])
rownames(protein.rpla.arraydat) = protein.rpla.f$Gene.name.c
colnames(protein.rpla.arraydat) = gsub('\\.','\\:',colnames(protein.rpla.f)[10:69])
protein.rpla.arraydat = protein.rpla.arraydat[match(unique(
 rownames(protein.rpla.arraydat)),rownames(protein.rpla.arraydat)),]

dna.arraydat = as.matrix(dnaf[,10:69])
rownames(dna.arraydat) = dnaf$Gene.name.c
colnames(dna.arraydat) = gsub('\\.','\\:',colnames(dnaf)[10:69])
dna.arraydat = dna.arraydat[match(unique(
 rownames(dna.arraydat)),rownames(dna.arraydat)),]
naCountByRow = apply(dna.arraydat,1,function(x) length(x[is.na(x)==T]))
dna.arraydat = dna.arraydat[naCountByRow<dim(dna.arraydat)[2],]

protein.ms.arraydat <- prep.protein.ms.assay(protein.ms.f)
deep.protein.ms.arraydat <- prep.protein.ms.assay(deep.protein.ms.f)
kinase.protein.ms.arraydat <- prep.protein.ms.assay(kinase.protein.ms.f)

# *********************************************** #
# format metadata associated with cell lines      #
# *********************************************** #
pdat = as.data.frame(cbind(pdataf$tissue.of.origin..a.,pdataf$age..a.,
 pdataf$sex..a.,pdataf$Contributor))
rownames(pdat) <- unlist(lapply(pdataf$Cell.Line.Name, 
 function(x) strsplit(x,' ')[[1]][1]))
colnames(pdat) <- c('tissue.of.origin','age','sex','contributor')

# ********************* #
# save data and pdata   #
# ********************* #
save(pdat,dna.arraydat,rna.arraydat,protein.rpla.arraydat,protein.ms.arraydat,
 deep.protein.ms.arraydat,file='nci60dat.rda')


