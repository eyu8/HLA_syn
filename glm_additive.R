library(readr)
library(data.table)

#arg 1 is FILE arg2 is REGION
args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1]
REGION <- args[2]

#load data depending if its UKB or not
if(FILE != "ukbb"){

    covar <- as.data.frame(fread(paste0("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/txt_data/", FILE, "/", FILE, "_covar.txt")))

    A <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-A_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    B <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-B_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    C <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-C_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DPB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DPB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DQA1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQA1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DQB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DRB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DRB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
} else if(FILE == "ukbb" && REGION == "PD"){

    PD <- as.data.frame(fread("~/runs/eyu8/data/ukbb_dagher/ukbb_PD_covar.txt"))
    control_PD <- as.data.frame(fread("~/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_control_PD_covar.txt"))

    PD$phenotype <- 2
    control_PD$phenotype <- 1

    covar <- rbind(PD, control_PD)

    A <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-A_ukbb_imp_HLA_Euro.csv"))
    B <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-B_ukbb_imp_HLA_Euro.csv"))
    C <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-C_ukbb_imp_HLA_Euro.csv"))
    DPB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DPB1_ukbb_imp_HLA_Euro.csv"))
    DQA1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQA1_ukbb_imp_HLA_Euro.csv"))
    DQB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQB1_ukbb_imp_HLA_Euro.csv"))
    DRB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DRB1_ukbb_imp_HLA_Euro.csv"))
    DRB3 <- as.data.frame(fread("ukbb_imp_HLA_DRB3/HLA-DRB3_ukbb_imp_HLA_DRB3.csv"))
    DRB4 <- as.data.frame(fread("ukbb_imp_HLA_DRB4/HLA-DRB4_ukbb_imp_HLA_DRB4.csv"))
    DRB5 <- as.data.frame(fread("ukbb_imp_HLA_DRB5/HLA-DRB5_ukbb_imp_HLA_DRB5.csv"))


} else if(FILE == "ukbb" && REGION == "Proxy"){

    Proxy <- as.data.frame(fread("~/runs/eyu8/data/ukbb_dagher/ukbb_proxy_covar.txt"))
    control_Proxy <- as.data.frame(fread("~/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_control_proxy_covar.txt"))

    Proxy$phenotype <- 2
    control_Proxy$phenotype <- 1

    covar <- rbind(Proxy, control_Proxy)

    A <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-A_ukbb_imp_HLA_Euro.csv"))
    B <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-B_ukbb_imp_HLA_Euro.csv"))
    C <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-C_ukbb_imp_HLA_Euro.csv"))
    DPB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DPB1_ukbb_imp_HLA_Euro.csv"))
    DQA1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQA1_ukbb_imp_HLA_Euro.csv"))
    DQB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQB1_ukbb_imp_HLA_Euro.csv"))
    DRB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DRB1_ukbb_imp_HLA_Euro.csv"))
    DRB3 <- as.data.frame(fread("ukbb_imp_HLA_DRB3/HLA-DRB3_ukbb_imp_HLA_DRB3.csv"))
    DRB4 <- as.data.frame(fread("ukbb_imp_HLA_DRB4/HLA-DRB4_ukbb_imp_HLA_DRB4.csv"))
    DRB5 <- as.data.frame(fread("ukbb_imp_HLA_DRB5/HLA-DRB5_ukbb_imp_HLA_DRB5.csv"))
}

#Merge allele with covariates
HLA <- list(A=merge(covar,A),B=merge(covar,B),C=merge(covar,C),
    DPB1=merge(covar,DPB1),DQA1=merge(covar,DQA1),DQB1=merge(covar,DQB1),DRB1=merge(covar,DRB1))

#Exclude samples with 2 alleles with posterior probability less than 0.5
poorSamples <- lapply(HLA,function(x) x[x$prob<0.5,]$sample.id)
poorSamples <- unlist(poorSamples, use.names=FALSE)
poorSamples <- poorSamples[duplicated(poorSamples)]
filtered_HLA <- lapply(HLA,function(x) x[!(x$sample.id %in% poorSamples),])

for(i in 7){
    message(paste0("Analysing HLA-",names(HLA)[i]))
    #Exclude allele where posterior probability less than 0.5
    filtered_HLA_gene <- filtered_HLA[[i]][filtered_HLA[[i]]$prob > 0.5,]
    #Extract alleles
    HLA_allele1 <- filtered_HLA_gene[, c("sample.id","allele1")]
    HLA_allele2 <- filtered_HLA_gene[, c("sample.id","allele2")]
    names(HLA_allele1) <- c("sample.id","allele")
    names(HLA_allele2) <- c("sample.id","allele")
    #Marge allele into 1 column
    HLA_allele <- rbind(HLA_allele1, HLA_allele2)
    HLA_count <- as.data.frame.matrix(table(HLA_allele))
    HLA_allele_name <- names(HLA_count)
    HLA_count$sample.id <- rownames(HLA_count)
    #Change homozygote to heterozygote for an additive model
    HLA_meta_data <- filtered_HLA_gene[, !(names(filtered_HLA_gene) %in% c("allele1", "allele2", "prob", "matching"))]
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
        if(FILE != "ukbb"){

            f <- formula(paste("phenotype ~ h + age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"))

        } else if(FILE == "ukbb"){

            f <- formula(paste("phenotype ~ h + Townsend + AgeAtRecruit + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"))

        }

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

        if(FILE == "ukbb" && REGION == "Proxy"){

            var_list <- append(var_list[1], var_list)
            var_list[[2]] <- var_list[[2]] * c(2,2,1)
            names(var_list[[2]]) <- c("h_b_adjusted", "h_StdErr_adjusted", "h_p_adjusted")
        }
        allele <- paste0(" ",allele)
        return(cbind(allele, HLA_ntotal, HLA_ncase, HLA_ncontrol, HLA_maf, HLA_maf_freq_case, HLA_maf_freq_control, Reduce(cbind, var_list)))})
    result <- Reduce(rbind, HLA_table)
    write_delim(result,paste0("results_",FILE,"/HLA-",names(HLA)[i],"_",FILE,"_",REGION,".csv"), delim = ",", na = "NA", append = FALSE, col_names = TRUE, escape = "double")
}
