
# This script assembles the assays on DNA, RNA and protein levels 
# on the NCI-60 panel cell lines into a multiAssayExperiment object. 
# Proteomics data cpome into two distict assays which vary the number of
# samples profiled and by the depth reached by the proteomic profiles. 

# The sample map objects store severla information not only on
# mapping samples from primary to assay-specic names bit also
# informatio on assays themselves.

require('MultiAssayExperiment')
load('nci60dat.rda')


rnadat <- Biobase::ExpressionSet(assayData=rna.arraydat)
proteinMSdat <- Biobase::ExpressionSet(assayData=protein.ms.arraydat)
deepProteinMSdat <- Biobase::ExpressionSet(assayData=deep.protein.ms.arraydat)
dnadat <- Biobase::ExpressionSet(assayData=dna.arraydat)

ExpList <- list(
 "expr1" = dnadat,
 "expr2" = rnadat,
 "expr3" = proteinMSdat,
 "expr4" = deepProteinMSdat
)
names(ExpList) <- tolower(names(ExpList))

# sample map in the DNA assay
dnamap = as.data.frame(cbind(rownames(pdat),colnames(dna.arraydat),
 rep("expr1",dim(pdat)[1]),rep("DNA:aCGH Agilent 44K",dim(pdat)[1]),
 rep("ratio of sample vs control of DNA copy number, normal female DNA 46,XX genomic DNA was obtained from Promega (Madison, WI)",dim(pdat)[1]),
 rep("Agilent Feature Extraction done with software version 8.1 with default settings for CGH arrays",dim(pdat)[1]),
 rep("http://discover.nci.nih.gov/cellminer/loadDownload.do",dim(pdat)[1])))
colnames(dnamap) = c('primary','assay','assayname','assaydescription',
 'rawdata','processingmethod',
 'source')
dnamap = dnamap[!is.na(dnamap$assay), ]

# sample map in the RNA assay
rnamap = as.data.frame(cbind(rownames(pdat),colnames(rna.arraydat),
 rep("expr2",dim(pdat)[1]),rep("RNA:Affy HG-U133 Plus 2.0",dim(pdat)[1]),
 rep("CEL file",dim(pdat)[1]),
 rep("RMA",dim(pdat)[1]),
 rep("http://discover.nci.nih.gov/cellminer/loadDownload.do",dim(pdat)[1])))
colnames(rnamap) = c('primary','assay','assayname','assaydescription',
 'rawdata','processingmethod','source')
rnamap = rnamap[!is.na(rnamap$assay), ]

# sample map in the MS-based proteomics assay
phenoCellLines = gsub('_','',unlist(lapply(rownames(pdat), 
 function(x) strsplit(x,":")[[1]][2])))
assayCellLines = unlist(lapply(colnames(protein.ms.arraydat), 
 function(x) strsplit(x,"_")[[1]][3]))
tmp = as.data.frame(cbind(rownames(pdat)[match(assayCellLines,phenoCellLines)], 
 colnames(protein.ms.arraydat),rep("expr3",length(assayCellLines)),
 rep("PROTEIN: mass spectrometry",length(assayCellLines)),
 rep("MS spectra",length(assayCellLines)),
 rep("MS spectra searched with Andromeda against IPI human database (maximum false discovery rate of 1%; protein abundance estimated based on summed peptide intensities)",length(assayCellLines)),
 rep("PMID:23933261 - Global proteome analysis of the NCI-60 cell line panel",
 length(assayCellLines))))
tmp = tmp[is.na(match(assayCellLines,phenoCellLines))==F,]
colnames(tmp) = c('primary','assay','assayname','assaydescription',
 'rawdata','processingmethod','source')
proteinMSmap = tmp

# sample map in the MS-based deep proteomics assay
phenoCellLines = gsub('_','',unlist(lapply(rownames(pdat),
 function(x) strsplit(x,":")[[1]][2])))
assayCellLines = unlist(lapply(colnames(deep.protein.ms.arraydat),
 function(x) strsplit(x,"_")[[1]][3]))
tmp = as.data.frame(cbind(rownames(pdat)[match(assayCellLines,phenoCellLines)],
 colnames(deep.protein.ms.arraydat),
 rep("expr4",length(assayCellLines)),rep("PROTEIN: deep mass spectrometry",
 length(assayCellLines)),
 rep("MS spectra",length(assayCellLines)),
 rep("MS spectra searched with Andromeda against IPI human database (maximum false discovery rate of 1%; protein abundance estimated based on summed peptide intensities)",length(assayCellLines)),
 rep("PMID:23933261 - Global proteome analysis of the NCI-60 cell line panel",
 length(assayCellLines))))
tmp = tmp[is.na(match(assayCellLines,phenoCellLines))==F,]
colnames(tmp) = c('primary','assay','assayname','assaydescription',
 'rawdata','processingmethod','source')
deepProteinMSmap = tmp


listmap <- list(dnamap,rnamap,proteinMSmap,deepProteinMSmap)
names(listmap) <- tolower(c("expr1","expr2","expr3","expr4"))
listmap

dfmap <- listToMap(listmap)
dfmap

names(ExpList) %in% as.character(dfmap[["assayname"]])

# remove unmapped sampes from assays
ExpList2 = list()
for(i in 1 : length(ExpList)){
 ExpList2[[i]] <- ExpList[[i]][,colnames(ExpList[[i]]) %in% listmap[[i]]$assay]
}
names(ExpList2) <- names(ExpList)

# check again
res <- lapply(1:length(listmap), function(i){
    has <- colnames(ExpList2[[i]]) %in% listmap[[i]]$assay
    colnames(ExpList2[[i]])[!has]
})
res


myMultiAssay <- MultiAssayExperiment(Elist=ExpList2, 
 pData=pdat, sampleMap=dfmap)


