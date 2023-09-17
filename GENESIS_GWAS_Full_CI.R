
##################################
# Load libraries
library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(SeqArray)
library(BiocParallel)
library(GENESIS)

# Convert PLINK data to .gds format
fileroot <- "All_CI"
snpgdsBED2GDS(bed.fn=paste0(fileroot,".bed"),
              bim.fn=paste0(fileroot,".bim"),
              fam.fn=paste0(fileroot,".fam"),
              out.gdsfn="All_CI.gds",family=TRUE)

# Open .gds file in SNPRelate
genofile <- snpgdsOpen("All_CI.gds")

samp.id <- read.gdsn(index.gdsn(genofile,"sample.id"))
fam.id <- read.gdsn(index.gdsn(genofile,"sample.annot/family"))

# Create and save kinship matrix 
ibd.robust <- snpgdsIBDKING(genofile,maf=0.01,missing.rate=0.02,
                            type="KING-robust",family.id=fam.id)
kinship <- ibd.robust$kinship
save(kinship,file="All_CI_kinship.Robj")
# Order of individuals in kinship matrix
save(samp.id,file="All_CI_kinship.sample.id.Robj")
snpgdsClose(genofile)

# Calculate PCs using GENESIS/PCAiR
kinship <- get(load("All_CI_kinship.Robj"))
kinship <- as.data.frame(kinship)
samp.id <- get(load("All_CI_kinship.sample.id.Robj"))
row.names(kinship) <- samp.id
colnames(kinship) <- samp.id 
genofile <- GdsGenotypeReader("All_CI.gds")
genodata <- GenotypeData(genofile)

# Get kinship matrix from KING
#file.king = system.file("All_CI_king_table.kin0")
#KINGmat <- kingToMatrix(
#  c(system.file("extdata", "All_CI_king_table.kin0", package="GENESIS"),
#    system.file("extdata", "All_CI.kin", package="GENESIS")),
 # estimator = "Kinship")
#KINGmat[1:5,1:5]

file.king = system.file("extdata", "All_CI_king_table.kin0", package = "GENESIS")
KINGmat <- kingToMatrix("All_CI_king.kin", estimator = "Kinship")

# Generate 20 PCs - save for later
mypcair <- pcair(genodata, kinobj = KINGmat, divobj = KINGmat)
pcs <- mypcair$vectors
pcsmat <- as.data.frame(cbind(row.names(pcs),pcs))
names(pcsmat) <- c("ID",paste0("PC",1:20))
write.table(pcsmat, "All_CI_PCAiR.txt",row.names=F,quote=F,sep="\t")

# Using PCRelate to generate kinship coefficients, conditional on 
# population structure. Unless specified, the file generated will 
# be called "tmp_pcrelate.gds"; use gds.prefix to specify file root name.
# When you choose write.to.gds=TRUE, then the file is saved as above. 
# If you say write.to.gds=FALSE then you get a pcrelate object.
genodata_block = GenotypeBlockIterator(genodata)
mypcrel <- pcrelate(genodata_block, pcs=pcs[,1:2], training.set=mypcair$unrels, BPPARAM = BiocParallel::SerialParam())
kinmat <- mypcrel$kinship
write.table(kinmat,"All_CI_PCRelate.txt.gz",row.names=F, quote=F, sep="\t")

#close(genodata) #Check if closed will fix later issue

# Save for later, so don't need to redo this part for other association tests using this dataset.
save(mypcrel, file=paste0("All_CI_mypcrel.RData"))


#==============================================================================
#Run GENESIS assocTestMM()
#
# mypcair contains PCs from a PC-AiR analysis.
# mypcrel contains Kinship Estimates (conditional on population structure using PC-AiR analysis) from a PC-Relate analysis.
# mypheno is a dataframe of Phenotype values.
#  - Make sure the order of individual is the same as in your genodata for pcair() and pcrelate().
#-------------

# phenotype file
pfile = "PV_Data_CI.csv"
mypheno = read.csv(pfile, header=TRUE, na.strings = "NA", stringsAsFactors = FALSE, as.is=T)
ind_order = read.table("All_CI.fam", head=F)
mypheno = mypheno[match(ind_order$V2, mypheno$IID),]
mypheno$sex[mypheno$sex == 1] = "M"
mypheno$sex[mypheno$sex == 2] = "F"

head(mypheno)

main_pheno = "dx"

# Get PC values
pcs = read.table("All_CI_PCAiR.txt", na.strings = "", header=T, stringsAsFactors = FALSE, as.is=T)

# Get GRM from mypcrel
myGRM = pcrelateToMatrix(mypcrel)

# Make analysis data.frame
mydat = data.frame(scanID = mypcair$sample.id, pc1 = pcs$PC1, pc2 = pcs$PC2,
                   pheno = mypheno$dx, age = mypheno$Age, sex = mypheno$sex)
head(mydat)

# Make ScanAnnotationDataFrame
scanAnnot = ScanAnnotationDataFrame(mydat)

# Fit the null mixed model
nullmod = fitNullModel(scanAnnot, outcome = "pheno", covars = c("age","sex","pc1","pc2"), 
                    cov.mat = myGRM, family = "binomial") #binary trait
#                    covMatList = myGRM) #quantitative trait

#genodata <- GenotypeData(HapMap_geno, scanAnnot = scanAnnot)
#genodata
genoIterator <- GenotypeBlockIterator(genodata, snpBlock=5000)

# Run single-SNP association test
assoc = assocTestSingle(genoIterator, null.model = nullmod, BPPARAM = BiocParallel::SerialParam())

save(assoc, file=paste0("All_CI_GENESIS_", main_pheno, ".RData"))

write.csv(assoc, file = "All_CI_GENESIS_GWAS.csv", quote=F, row.names=F)
