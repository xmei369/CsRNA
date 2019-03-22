setwd('C:/Users/Xin Mei/Desktop/gadProperty')

CountTable = read.csv("featureCounts_matrix.csv", header=TRUE, row.names=1)
head(CountTable)

'
Design = data.frame(
  row.names = colnames( CountTable ),
  condition = c( "early_stage_lateral_bud", "apical_bud", "lateral_bud",
                 "the_first_leaf", "the_second_leaf", "mature_leaf", "old_leaf",
                 "one_and_a_bud", "two_and_a_bud", "root", "stem", "flower", "seed") )
'
Design = data.frame(
  row.names = colnames( CountTable ),
  condition = c( "ELB", "ApiB", "LatB", "1stL", "2ndL", "matL", "oldL",
                 "1B1L", "1B2L", "Ro", "St", "Fl", "Se") )
Design

condition = Design$condition
colnames(CountTable) <- condition

# CsGAD expr level:
'
TEA021063.1: GAD1(short) 208960 
TEA021095.1: GAD1 209300 
TEA024088.1: GAD2 168685 
'
gad <- CountTable[rownames(CountTable) %in% c("TEA021095.1", "TEA024088.1"),]
gad
barplot(as.matrix(gad[1,]))
dotplot(as.matrix(gad[2,]))

# CsCaM expr level:
"grep 'calmodulin;' C.sinensis.gene.kegg.annotation.xls"
CaM <- CountTable[rownames(CountTable) %in% c("TEA021147.1","TEA018823.1", "TEA031520.1"),]
CaM
barplot(as.matrix(CaM[2,]))


library(DESeq)
cds = newCountDataSet( CountTable, condition )

cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )

# CsGAD normalized expr level:
gad <- counts( cds, normalized=TRUE )[rownames(CountTable) %in% c("TEA021095.1", "TEA024088.1"),]
gad
barplot(gad[1,])
barplot(gad[2,])
write.csv(gad, "gadNormExpr.csv")

'
cds = estimateDispersions( cds, method="blind", sharingMode="fit-only" )
str( fitInfo(cds) )
plotDispEsts( cds )
head( fData(cds) )

res = nbinomTest( cds, "two_and_a_bud", "stem" )
head(res)
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")

res[res$id %in% c("TEA021063.1","TEA021095.1","TEA024088.1"),]
'

#resSig = res[ res$padj < 0.1, ]
#head( resSig[ order(resSig$pval), ] )
#head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )
#head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )
#write.csv( res, file="My Pasilla Analysis Result Table.csv" )
