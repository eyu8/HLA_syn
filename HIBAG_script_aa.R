#module load gcc/9.3 r-bundle-bioconductor/3.12

library(data.table)

library(HIBAG)
library(readr)

#arg 1 is FILE arg2 is REGION
args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1]
REGION <- args[2]

covar <- as.data.frame(fread(paste0("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/txt_data/", FILE, "/", FILE, "_covar.txt")))

A <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-A_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
B <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-B_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
C <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-C_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
DPB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DPB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
DQA1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQA1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
DQB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
DRB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DRB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)

HLA <- list(A=merge(covar,A),B=merge(covar,B),C=merge(covar,C),
            DPB1=merge(covar,DPB1),DQA1=merge(covar,DQA1),DQB1=merge(covar,DQB1),DRB1=merge(covar,DRB1))

for(i in 1:7){
    message(paste0("Analysing HLA-",names(HLA)[i]))
    hla_group <- HLA[[i]]
    hla_group$phenotype[hla_group$phenotype == -9] <- NA
    hla_group$phenotype <- as.factor(hla_group$phenotype - 1)
    hla <- hlaAllele(hla_group$"sample.id", H1=hla_group$allele1, H2=hla_group$allele2, locus=names(HLA[i]), assembly="hg19", prob=hla_group$prob)
    name <- names(HLA[i])
    aa <- hlaConvSequence(hla, code="P.code.merge")
    hla_aa <- aa$value
    hla_aa_table <- summary(aa)

    #Increment position
    increment <- hla_aa_table[1,"Pos"] - 1
    hla_aa_table[,"Pos"] <- hla_aa_table[,"Pos"] - increment

    hla_pos <- hla_aa_table[,"Pos"]
    HLA_allele1 <- HLA[, c("sample.id","allele1")]
    HLA_allele2 <- HLA[, c("sample.id","allele2")]
    names(HLA_allele1) <- c("ID","allele")
    names(HLA_allele2) <- c("ID","allele")
    HLA_allele <- rbind(HLA_allele1, HLA_allele2)

    hla_aa_result <- lapply(hla_pos, function(pos){
            HLA_aa_code <- HLA_allele
            HLA_aa_code$allele <- substr(HLA_aa_code$allele, pos, pos)
            HLA_count <- as.data.frame.matrix(table(HLA_aa_code))
            names(HLA_count)[names(HLA_count) == "-"] <- "Ref"
            names(HLA_count)[names(HLA_count) == "*"] <- "Amb"
            names(HLA_count)[names(HLA_count) == "."] <- "CNV"
            HLA_allele_name <- names(HLA_count)
            HLA_count$"sample.id" <- rownames(HLA_count)
            if(any(HLA_count$All == 0)){
                HLA_count[HLA_count$All == 0,]$All <- 1
            }
            if(any(HLA_count$All == 2)){
                HLA_count[HLA_count$All == 2,]$All <- 0
            }
            HLA_meta_data <- hla_group[, !(names(hla_group) %in% c("allele1", "allele2", "prob", "matching"))]
            HLA_glm <- merge(HLA_meta_data, HLA_count)
            HLA_table <- lapply(HLA_allele_name, function(allele){
                HLA_glm$sex <- as.factor(HLA_glm$sex)
                HLA_glm[,allele] <- as.numeric(HLA_glm[,allele])

                if(length(levels(as.factor(HLA_glm[,allele]))) == 1){
                        return(NULL)
                }
                f <- formula(paste("phenotype ~ ",allele," + age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"))
                HLA_result <- glm(f,family = binomial, data = HLA_glm)
                HLA_carrier_freq <-  table(HLA_glm[, c(allele, "phenotype")])
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
                CI <- signif(as.data.frame(confint.default(HLA_result)),3)
                var <- signif(as.data.frame(coef(summary(HLA_result))),3)
            if(!grepl(allele,row.names(var)[2])){
                return(NULL)
            }
                var_list <- lapply(2:nrow(var), function(i){
                    summary_stats <- cbind(var[i, 1, drop = FALSE], var[i, 2, drop = FALSE], var[i, 4, drop = FALSE])
                    names(summary_stats) <- c("b", "StdErr", "p")
                    if(i == 2){
                        row.names(summary_stats) <- "aa"
                        }
                    names(summary_stats) <- paste(row.names(summary_stats), names(summary_stats), sep = "_")
                    return(summary_stats)
                    })
                return(cbind(Pos = paste0("Pos",pos + increment,allele), HLA_ntotal, HLA_ncase, HLA_ncontrol, HLA_maf, HLA_maf_freq_case, HLA_maf_freq_control, Reduce(cbind, var_list)))})
            result <- Reduce(rbind, HLA_table)
            return(result)
            })

    hla_result <- Reduce(rbind, hla_aa_result)
    write_delim(hla_result,paste0("results_",FILE,"/HLA-",names(HLA)[i],"_",FILE,"_",REGION,"_aa.csv"), delim = ",", na = "NA", append = FALSE, col_names = TRUE, escape = "double")
}
