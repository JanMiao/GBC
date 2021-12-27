# GBC
Genomic Breed Composition analysis and breed identification

## 联系方式
miaojian@zju.edu.cn

## Dependency
**R package**
1. data.table
2. ggplot2
3. glmnet (for lasso analysis)
4. quadprog (for linear model with constrained regression coefficients)

**Plink** (>=1.9)

## Install 
install.packages("Path_to/GBC_1.0_R_x86_64-pc-linux-gnu.tar.gz", repos = NULL)

## Quickly start
快速生成参考群，并基于该参考群进行品种溯源分析和品种鉴定分析，使用如下三个函数：
- `Makeref` 生成参考群
- `GBCpred` GBC评估
- `BIpred` 品种鉴定
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
基于plink 二进制（bed）文件和样本对应的品种信息，对样本进行outliers/relatives deletion。去除异常样本（错标/杂种）并对高亲缘关系的样本进行去冗余。  

**Usage**     
`Makeref(plink_prefix, breed_file, quant_thres,
        IF_RMrelates, relates_thres, plink_dir, REF_prefix="REF") `

**Return**   
当前目录下生成**REF_frequency.rds** 和 **REF_genotype.rds**两个文件。**REF_frequency.rds**是每个品种在每个位点的等位基因频率文件。**REF_genotype.rds**是一个list，包含两个元素。其一是质控后样本的genotype文件，已经过plink `recodeA`（编码为0,1,2）；其二是质控后样本ID及其对应的品种信息。

**Parameter**   
`plink_prefix`：plink 二进制文件前缀。其中fam文件第二列必须是品种ID，且没有重复。 
`breed_file`：品种信息文件。没有表头，第一列是品种ID，第二列是对应是品种名（英文缩写，课有英文字母、数字、下划线组成）。 
`quant_thres`：保留多少比例的样本。程序根据genotype likelihood排序所有样本。如`quant_thres=0.9`， 则保留genotype likelihood排序前90%的样本。其余10%的样本视为outliers。 
`IF_RMrelates`: 是否需要剔除亲缘关系系数过高的样本。
`relates_thres`: 亲缘关系系数阈值。仅当`IF_RMrelates=True`时生效。 
`plink_dir`: plink 二进制文件所在目录。
`REF_prefix`: 生成的两个rds文件的前缀。

**Example**       
`Makeref(plink_prefix="ref", breed_file="breedInfo.txt", quant_thres = 0.9,
        IF_RMrelates=TRUE, relates_thres = 0.8, plink_dir=plink_dir, REF_prefix="REF")`

### `GBCpred`
**Description**   
依据`Makeref`生成的**REF_frequency.rds**，对test data 进行GBC analysis。  

**Usage**     
`GBCpred(RDS, test_prefix, method, plink_dir, breedused=NULL, nmarker=NULL, BINnum=500, S=0.01, testMode=FALSE)`  

**Return**   
每个样本的GBC分析结果。  

**Parameter**   
`RDS`: 等位基因频率RDS文件名。 
`test_prefix`：代评估plink 二进制文件前缀。其中fam文件第二列必须是品种ID，且没有重复。  
`method`：*lm, lasso, lmC, similar* 中的一个。 
`plink_dir`：plink 二进制文件所在目录。 
`breedused`：GBC分析所用的品种。默认`NULL`代表使用所有品种。 
`nmarkers`：breed information marker的数目。默认`NULL`代表使用所有marker。 
`BINnum`：单个滑块包含的marker数目。 仅当`method=similar`时生效。 
`S`：LASSO的lambda参数值，lambda越大，LASSO对回归系数的惩罚越强。  
`testMode`：如设置为`TRUE`，会再当前目录下生成**test_prefixGBCdata.rds**文件，方便对相同的test data进行多次尝试不同模型或参数。

**Example**    
`gbc1 = GBCpred(RDS="REF_frequency.rds",
               test_prefix="test", nmarkers=NULL, testMode = TRUE,
               method="lm", plink_dir=plink_dir)`


### `BIpred`
**Description**   
依据`Makeref`生成的**REF_frequency.rds**，对test data 进行品种鉴定分析。  

**Usage**   
`BIpred(RDS, test_prefix, plink_dir, nmarkers=NULL, breedused=NULL, testMode=FALSE)`  

**Return**    
每个样本的品种鉴定分析结果。   

**Parameter**   
`RDS`: 等位基因频率RDS文件名。 
`test_prefix`：代评估plink 二进制文件前缀。其中fam文件第二列必须是品种ID，且没有重复。 
`plink_dir`：plink 二进制文件所在目录。 
`breedused`：GBC分析所用的品种。默认`NULL`代表使用所有品种。 
`nmarkers`：breed information marker的数目。默认`NULL`代表使用所有marker。   
`testMode`：如设置为`TRUE`，会再当前目录下生成**test_prefixGBCdata.rds**文件，方便对相同的test data进行多次尝试不同模型或参数。

**Example**     
`bi = BIpred(RDS="REF_frequency.rds", test_prefix="test", nmarkers=20, plink_dir=plink_dir)`  


### `GBCplot`
**Description**   
将GBCpred生成的GBC结果进行可视化。

**Usage**   
`GBCplot(GBCres, FontSize)`  

**Return**   
GBC可视化结果的ggplot2对象。 

**Parameter**   
`GBCres`: GBCpred得到的GBC结果。 
`FontSize`：样本ID的字体大小。

**Example**     
`p = GBCplot(GBCres=GBCmatrix, FontSize=2)`  
`p`


### `PCAplot`
**Description**   
绘制PCA二维图。

**Usage**   
`PCAplot(plink_prefix, sampleID, breed_file, PDFprefix, plink_dir)`  

**Return**   
在当前目录下生成PCA二维图。 

**Parameter**   
`plink_prefix`：plink 二进制文件前缀。其中fam文件第二列必须是品种ID，且没有重复。 `sampleID`：需要绘图的样本ID。

`breed_file`：品种信息文件。没有表头，第一列是品种ID，第二列是对应是品种名（英文缩写，课有英文字母、数字、下划线组成）。 
`PDFprefix`：生成的pdf文件的前缀。 
`plink_dir`: plink二进制文件所在目录。

**Example**   
`PCAplot(plink_prefix="Chip", breed_file="breedChip.txt", PDFprefix="QC",
        plink_dir="~/software/")`

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

