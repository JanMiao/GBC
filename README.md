# GBC
Genomic Breed Composition analysis and breed identification

## contact information
miaojian6363@163.com

## Dependency
**R package**
1. data.table
2. ggplot2
3. glmnet (for lasso analysis, **optional**)
4. quadprog (for linear model with constrained regression coefficients, **optional**)

**Plink** (>=1.9)

## Install 
install.packages("Path_to/GBC_1.0_R_x86_64-pc-linux-gnu.tar.gz", repos = NULL)

## Quickly start
Performing Genomic Breed Composition analysis and breeding identification using reference dataset and generating reference dataset using user-provided data. We also provide a preprrocessed reference dataset of some pig breeds. This R package has below three main functions:
- `Makeref` construction of reference data
- `GBCpred` GBC analysis
- `BIpred` breed identification
```
### copy example files to your current work directory
example_files=c(
  system.file("extdata", "ref.bed", package = "GBC"),
  system.file("extdata", "ref.bim", package = "GBC"),
  system.file("extdata", "ref.fam", package = "GBC"),
  system.file("extdata", "test.bed", package = "GBC"),
  system.file("extdata", "test.bim", package = "GBC"),
  system.file("extdata", "test.fam", package = "GBC"),
  system.file("extdata", "breedInfo.txt", package = "GBC")
)
file.copy(example_files, ".")

### specify your directory of plink executable file
plink_dir = "~/software/"

### 1. generate your reference data
Makeref(plink_prefix="ref", breed_file="breedInfo.txt", quant_thres = 0.9,
        IF_RMrelates=TRUE, relates_thres = 0.8,
        plink_dir=plink_dir)

### 2. Predict GBC of test data Using reference data
# lm
gbc1 = GBCpred(RDS="REF_frequency.rds",
               test_prefix="test", nmarkers=NULL, testMode = T,
               method="lm", plink_dir=plink_dir)
head(gbc1)
p = GBCplot(gbc1, FontSize=4)

#lasso
require(glmnet)
gbc2 = GBCpred(RDS="REF_frequency.rds",
               test_prefix="test", nmarkers=NULL, testMode = T,
               method="lasso", plink_dir=plink_dir)

### 3. Predict Breed of test data Using reference data
bi = BIpred(RDS="REF_frequency.rds",
            test_prefix="test", nmarkers=20,
            plink_dir=plink_dir)
```

## Function description
### `Makeref`
**Description**   
Remove outliers and relatives in user-provided plink files and then generate two rds file for downstram GBC/BI analysis.

**Usage**     
`Makeref(plink_prefix, breed_file, quant_thres,
        IF_RMrelates, relates_thres, plink_dir, REF_prefix="REF") `

**Return**   
**REF_frequency.rds** and **REF_genotype.rds** are generate in current working dictionary。**REF_frequency.rds** is the allele frequency file。**REF_genotype.rds**is a list including two elements, one of which is the plink `recodeA` genotypes(coding as 0, 1, 2) and the other is the sample IDs and the corresponding breed labels.

**Parameter**   
`plink_prefix`：The prefix of plink bed file. The second column of fam file is the sample ID.
`breed_file`：File contains breed information, no header，the first column is the sample ID, the second column is the breed label.
`quant_thres`：Propotion of sample should be treated as pure. if `quant_thres=0.9`, then keep 90% pure samples according the order of genotype likelihood. 
`IF_RMrelates`: If relatives should be removed.
`relates_thres`: The threshold of coefficient of relatedness. This parameter works only when `IF_RMrelates=True`.
`plink_dir`: The dictionary of plink executable file.
`REF_prefix`: The prefix of two rds files.

**Example**       
`Makeref(plink_prefix="ref", breed_file="breedInfo.txt", quant_thres = 0.9,
        IF_RMrelates=TRUE, relates_thres = 0.8, plink_dir=plink_dir, REF_prefix="REF")`

### `GBCpred`
**Description**   
Perform GBC analysis, according to **REF_frequency.rds** generated by `Makeref`

**Usage**     
`GBCpred(RDS, test_prefix, method, plink_dir, breedused=NULL, nmarker=NULL, BINnum=500, S=0.01, testMode=FALSE)`  

**Return**   
GBC of all testing samples. 

**Parameter**   
`RDS`: The allele frequence file.
`test_prefix`：The prefix of plink bed file (testing samples). The second column of fam file is the sample ID. 
`method`：One from *lm, lasso, lmC, similar*  
`plink_dir`：The dictionary of plink executable file. 
`breedused`：A vector of subset breeds from breeds in reference data.
`nmarkers`：Number of breed information marker. `NULL` means using all marker. 
`BINnum`：The marker number in a single bin,  works only when `method=similar`. 
`S`：lambda value in LASSO.   
`testMode`：If `TRUE`, the intermediate file will be saved and the redundant steps will be skipped. Only used when testing data are consistent.

**Example**    
`gbc1 = GBCpred(RDS="REF_frequency.rds",
               test_prefix="test", nmarkers=NULL, testMode = TRUE,
               method="lm", plink_dir=plink_dir)`


### `BIpred`
**Description**   
Perform breed identification(BI) analysis, according to **REF_frequency.rds** generated by `Makeref`

**Usage**   
`BIpred(RDS, test_prefix, plink_dir, nmarkers=NULL, breedused=NULL, testMode=FALSE)`  

**Return**    
breed identification results. 

**Parameter**   
`RDS`: The allele frequence file.
`test_prefix`：The prefix of plink bed file (testing samples). The second column of fam file is the sample ID. 
`plink_dir`：The dictionary of plink executable file. 
`breedused`：A vector of subset breeds from breeds in reference data.
`nmarkers`：Number of breed information marker. `NULL` means using all marker.   
`testMode`：If `TRUE`, the intermediate file will be saved and the redundant steps will be skipped. Only used when testing data are consistent.

**Example**     
`bi = BIpred(RDS="REF_frequency.rds", test_prefix="test", nmarkers=20, plink_dir=plink_dir)`  


### `GBCplot`
**Description**   
plot GBC results.

**Usage**   
`GBCplot(GBCres, FontSize)`  

**Return**   
a object of ggplot2

**Parameter**   
`GBCres`: the results from GBCpred 
`FontSize`：font size.

**Example**     
`p = GBCplot(GBCres=GBCmatrix, FontSize=2)`  
`p`

## Detail functions
```
# preprocess: data preprocess
d1 = preprocess(plink_prefix="ref", breed_file="breedInfo.txt",
           plink_dir=plink_dir)

# genotypeLL: detect outliers for a population using genotype likelihood
p1_b = d1$breedInfo[d1$breedInfo=="Yorkshire"]
p1 = d1$plinkAmatrix[names(p1_b),]
lh = genotypeLL(p1)

# sampleRM: remove relatives based on genomic correlation value
g = makeG(p1)
sampleRM = names(RmRelates(g, 0.8))

# GetTestGenotype: obtain same markers among reference data and test data
sMarker = GetTestGenotype(freqMatrix=RefFreq, test_prefix="test",
                plink_dir=plink_dir)
names(sMarker)

# GBClm: GBC analysis using lm
gbc = GBClm(genotypeMatrix=sMarker$testGeno, AFMatrix=sMarker$RefFrqMatrix)
# choose possible breeds
gbc = GBClm(genotypeMatrix=sMarker$testGeno, AFMatrix=sMarker$RefFrqMatrix,
             breedused=c("DU","LR","LW"))

# GBClasso: GBC analysis using lasso
require(glmnet)
gbc = GBClasso(genotypeMatrix=sMarker$testGeno, AFMatrix=sMarker$RefFrqMatrix)

# GBClmConstrain: GBC analysis using lm with constrained parameters
require(quadprog)
gbc = GBClmConstrain(genotypeMatrix=sMarker$testGeno, AFMatrix=sMarker$RefFrqMatrix)

# GBCsimilarity: GBC analysis using genome similarity
gbc = GBCsimilarity(genotypeMatrix=sMarker$testGeno, AFMatrix=sMarker$RefFrqMatrix,
                    SNPnum=3)

# GBCadjust_set0: adjust GBC result
gbc = GBClm(genotypeMatrix=sMarker$testGeno, AFMatrix=sMarker$RefFrqMatrix)
gbc = GBCadjust_set0(gbc, Ndecimals = 3)

# LLBI: Breed Identification using genotype likelihood
bi = LLBI(genotypeMatrix=sMarker$testGeno, AFMatrix=sMarker$RefFrqMatrix)
```

