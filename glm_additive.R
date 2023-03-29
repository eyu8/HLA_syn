library(readr)
library(data.table)

#arg 1 is FILE arg2 is REGION
args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1]
REGION <- args[2]

#load data depending if its UKB or not

covar <- as.data.frame(fread(paste0("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/txt_data/", FILE, "/", FILE, "_covar.txt")))

A <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-A_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
B <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-B_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
C <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-C_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
DPB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DPB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
DQA1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQA1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
DQB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
DRB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DRB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)

#Merge allele with covariates
HLA <- list(A=merge(covar,A),B=merge(covar,B),C=merge(covar,C),
    DPB1=merge(covar,DPB1),DQA1=merge(covar,DQA1),DQB1=merge(covar,DQB1),DRB1=merge(covar,DRB1))

for(i in 1:7){
    message(paste0("Analysing HLA-",names(HLA)[i]))

    #Extract alleles
    HLA_allele1 <- HLA[, c("sample.id","allele1")]
    HLA_allele2 <- HLA[, c("sample.id","allele2")]
    names(HLA_allele1) <- c("sample.id","allele")
    names(HLA_allele2) <- c("sample.id","allele"
    )
    #Marge allele into 1 column
    HLA_allele <- rbind(HLA_allele1, HLA_allele2)
    HLA_count <- as.data.frame.matrix(table(HLA_allele))
    HLA_allele_name <- names(HLA_count)
    HLA_count$sample.id <- rownames(HLA_count)
    
    #Change homozygote to heterozygote for an additive model
    HLA_meta_data <- HLA[, !(names(HLA) %in% c("allele1", "allele2", "prob", "matching"))]
    #Merge back alleles
    HLA_glm <- merge(HLA_meta_data, HLA_count)
    #Set sex, phenotype as factor; exclude alleles not present in controls
    HLA_table <- lapply(HLA_allele_name, function(allele){
        HLA_glm$phenotype[HLA_glm$phenotype == -9] <- NA
        HLA_glm$phenotype <- as.factor(HLA_glm$phenotype - 1)
        HLA_glm$sex <- as.factor(HLA_glm$sex)
        names(HLA_glm)[names(HLA_glm) == allele] <- "h"
        HLA_glm$h <- as.numeric(HLA_glm$h)
        if(length(levels(as.factor(HLA_glm$h))) == 1 ){
            return(NULL)
        }
        if(length(levels(as.factor(HLA_glm$phenotype))) == 1){
            return(NULL)
        }
        f <- formula(paste("phenotype ~ h + age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"))

        HLA_result <- glm(f,family = binomial, data = HLA_glm)
        #Get frequency and count
        HLA_carrier_freq <-  table(HLA_glm[, c("h", "phenotype")])
        HLA_h0_p0 <- HLA_carrier_freq[1,1]
        HLA_h0_p1 <- HLA_carrier_freq[1,2]
        HLA_h1_p0 <- HLA_carrier_freq[2,1]
        HLA_h1_p1 <- HLA_carrier_freq[2,2]

        if(nrow(HLA_carrier_freq)==2){
        HLA_h2_p0 <- 0
        HLA_h2_p1 <- 0
        }else{
            HLA_h2_p0 <- HLA_carrier_freq[3,1]
            HLA_h2_p1 <- HLA_carrier_freq[3,2]
        }

        
        HLA_maf <- signif((HLA_h1_p0 + HLA_h1_p1 + 2*HLA_h2_p0+2*HLA_h2_p1)/(sum(HLA_carrier_freq)*2),3)
        HLA_maf_freq_case <- signif((HLA_h1_p1+2*HLA_h2_p1)/(2*HLA_h0_p1+2*HLA_h1_p1+2*HLA_h2_p1),3)
        HLA_maf_freq_control <- signif((HLA_h1_p0+2*HLA_h2_p0)/(2*HLA_h0_p0+2*HLA_h1_p0+2*HLA_h2_p0),3)
        HLA_ncase <- HLA_h0_p1 + HLA_h1_p1 + HLA_h2_p1
        HLA_ncontrol <- HLA_h0_p0 + HLA_h1_p0 + HLA_h2_p0
        HLA_ntotal <- sum(HLA_carrier_freq)
        #Get statistics (beta, se, pvalue)
        var <- signif(as.data.frame(coef(summary(HLA_result))),3)
        if(!grepl("h",row.names(var)[2])){
            return(NULL)
        }
        var_list <- lapply(2:nrow(var), function(i){
            summary_stats <- cbind(var[i, 1, drop = FALSE], var[i, 2, drop = FALSE], var[i, 4, drop = FALSE])
            names(summary_stats) <- c("b", "StdErr", "p")
            names(summary_stats) <- paste(row.names(summary_stats), names(summary_stats), sep = "_")
            return(summary_stats)
        })

        allele <- paste0(" ",allele)
        return(cbind(allele, HLA_ntotal, HLA_ncase, HLA_ncontrol, HLA_maf, HLA_maf_freq_case, HLA_maf_freq_control, Reduce(cbind, var_list)))})
    result <- Reduce(rbind, HLA_table)
    write_delim(result,paste0("results_",FILE,"/HLA-",names(HLA)[i],"_",FILE,"_",REGION,".csv"), delim = ",", na = "NA", append = FALSE, col_names = TRUE, escape = "double")
}
