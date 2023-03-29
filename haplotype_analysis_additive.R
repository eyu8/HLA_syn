#module load gcc/9.3 r-bundle-bioconductor/3.12

library(HIBAG)
library(readr)
library(data.table)
library(haplo.stats)

#arg 1 is FILE arg2 is REGION
args <- commandArgs(trailingOnly = TRUE)


#FILE <- "MCGILL"
#REGION <- "Euro"

FILE <- args[1]
REGION <- args[2]

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

    PD <- as.data.frame(fread("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_PD_covar.txt"))
    control_PD <- as.data.frame(fread("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_control_PD_covar.txt"))

    names(PD)[names(PD) == "AgeAtRecruit"] <- "age"
    names(control_PD)[names(control_PD) == "AgeAtRecruit"] <- "age"

    PD$phenotype <- 2
    control_PD$phenotype <- 1

    names(PD)[1] <- "sample.id"
    names(control_PD)[1] <- "sample.id"

    covar <- rbind(PD, control_PD)

    A <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-A_ukbb_imp_HLA_Euro.csv"))
    B <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-B_ukbb_imp_HLA_Euro.csv"))
    C <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-C_ukbb_imp_HLA_Euro.csv"))
    DPB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DPB1_ukbb_imp_HLA_Euro.csv"))
    DQA1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQA1_ukbb_imp_HLA_Euro.csv"))
    DQB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQB1_ukbb_imp_HLA_Euro.csv"))
    DRB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DRB1_ukbb_imp_HLA_Euro.csv"))

} else if(FILE == "ukbb" && REGION == "Proxy"){

    Proxy <- as.data.frame(fread("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_proxy_covar.txt"))
    control_Proxy <- as.data.frame(fread("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_control_proxy_covar.txt"))
    names(Proxy)[names(Proxy) == "AgeAtRecruit"] <- "age"
    names(control_Proxy)[names(control_Proxy) == "AgeAtRecruit"] <- "age"

    Proxy$phenotype <- 2
    control_Proxy$phenotype <- 1

    names(Proxy)[1] <- "sample.id"
    names(control_Proxy)[1] <- "sample.id"

    covar <- rbind(Proxy, control_Proxy)

    A <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-A_ukbb_imp_HLA_Euro.csv"))
    B <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-B_ukbb_imp_HLA_Euro.csv"))
    C <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-C_ukbb_imp_HLA_Euro.csv"))
    DPB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DPB1_ukbb_imp_HLA_Euro.csv"))
    DQA1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQA1_ukbb_imp_HLA_Euro.csv"))
    DQB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQB1_ukbb_imp_HLA_Euro.csv"))
    DRB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DRB1_ukbb_imp_HLA_Euro.csv"))
}

HLA <- list(A=merge(covar,A),B=merge(covar,B),C=merge(covar,C),
    DPB1=merge(covar,DPB1),DQA1=merge(covar,DQA1),DQB1=merge(covar,DQB1),DRB1=merge(covar,DRB1))

poorSamples <- lapply(HLA,function(x) x[x$prob<0.5,]$sample.id)
poorSamples <- unlist(poorSamples, use.names=FALSE)

filtered_HLA <- lapply(HLA,function(x) x[!(x$sample.id %in% poorSamples), c("sample.id","allele1","allele2")])


HLA.list <- lapply(names(filtered_HLA),function(gene){
        HLA_gene <- filtered_HLA[[gene]]
        names(HLA_gene) <- c("sample.id",paste(gene,c("a1","a2"),sep = "."))
        return(HLA_gene)
        })

HLA_merged <- merge(covar, Reduce(merge, HLA.list))
HLA_covar <- HLA_merged[, !(names(HLA_merged) %in% grep(".a",names(HLA_merged)))]



seed <- c(17, 53, 1, 40, 37, 0, 62, 56, 5, 52, 12, 1)
set.seed(seed)

haplotype_glm <- function(pattern){

    message(paste("Running", pattern))

    geno <- HLA_merged[,grep(pattern,names(HLA_merged))]
    label <- names(HLA)[grep(pattern,names(HLA))]

        geno.glm <- setupGeno(geno, miss.val=c(0,NA), locus.label=label)

        glm.data <- cbind(HLA_covar,data.frame(geno.glm))
        glm.data$phenotype[glm.data$phenotype == -9] <- NA
        glm.data$phenotype <- glm.data$phenotype - 1
        glm.data$sex <- factor(glm.data$sex)

        save.em <- haplo.em(geno=geno, locus.label=label)

        haplCalls <- data.table(sample.id=HLA_covar$sample.id[save.em$subj.id], hap.1=save.em$hap1code,  hap.2=save.em$hap2code, hap.prob=save.em$post)
        haplotypes <- apply(save.em$haplotype, 1, function(x) paste0(paste(names(x), x , sep = "_"), collapse = "_"))

        haplCalls$hap.1 <- haplotypes[haplCalls$hap.1]
        haplCalls$hap.2 <- haplotypes[haplCalls$hap.2]
        filtered_haplCalls <- subset(haplCalls, hap.prob > 0.2)

        em_allele1 <- filtered_haplCalls[, c("sample.id","hap.1")]
        em_allele2 <- filtered_haplCalls[, c("sample.id","hap.2")]
        names(em_allele1) <- c("sample.id","hap")
        names(em_allele2) <- c("sample.id","hap")
        em_allele <- rbind(em_allele1, em_allele2)
        em_count <- as.data.frame.matrix(table(em_allele))
        haploNames <- names(em_count)
        em_count$sample.id <- row.names(em_count)
        haplo_glm <- merge(glm.data, em_count)

        haplo_table <- lapply(haploNames, function(hap){
            if(length(levels(as.factor(haplo_glm[,hap]))) == 1 ){
            return(NULL)
            }
            names(haplo_glm)[names(haplo_glm) == hap] <- "haplo"
            if(FILE != "ukbb"){

                f <- formula("phenotype ~ haplo + age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")

            } else if(FILE == "ukbb"){

                f <- formula("phenotype ~ haplo + Townsend + age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")

            }

            HLA_result <- glm(f, family = binomial, data = haplo_glm)

            Haplo_ntotal <- nrow(haplo_glm)
            Haplo_ncase <- nrow(haplo_glm[haplo_glm$phenotype == 1,])
            Haplo_ncontrol <- nrow(haplo_glm[haplo_glm$phenotype == 0,])

            Haplo_h0_p0 <- nrow(subset(haplo_glm, phenotype==0 & haplo==0))
            Haplo_h0_p1 <- nrow(subset(haplo_glm, phenotype==1 & haplo==0))
            Haplo_h1_p0 <- nrow(subset(haplo_glm, phenotype==0 & haplo==1))
            Haplo_h1_p1 <- nrow(subset(haplo_glm, phenotype==1 & haplo==1))
            Haplo_h2_p0 <- nrow(subset(haplo_glm, phenotype==0 & haplo==2))
            Haplo_h2_p1 <- nrow(subset(haplo_glm, phenotype==1 & haplo==2))

            Haplo_maf <- signif((Haplo_h1_p0 + Haplo_h1_p1 + 2*Haplo_h2_p0+2*Haplo_h2_p1)/(nrow(haplo_glm))*2,3)
            Haplo_maf_freq_case <- signif((Haplo_h1_p1+2*Haplo_h2_p1)/(2*Haplo_h0_p1+2*Haplo_h1_p1+2*Haplo_h2_p1),3)
            Haplo_maf_freq_control <- signif((Haplo_h1_p0+2*Haplo_h2_p0)/(2*Haplo_h0_p0+2*Haplo_h1_p0+2*Haplo_h2_p0),3)
            var <- signif(as.data.frame(coef(summary(HLA_result))),3)
            if(!grepl("haplo",row.names(var)[2])){
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
                names(var_list[[2]]) <- c("b_adjusted", "StdErr_adjusted", "p_adjusted")
            }
            Haplotype <- hap
            return(cbind(Haplotype, Haplo_ntotal, Haplo_ncase, Haplo_ncontrol, Haplo_maf, Haplo_maf_freq_case, Haplo_maf_freq_control, Reduce(cbind, var_list)))})



        result <- Reduce(rbind, haplo_table)
        filename <- paste(label, collapse = "_")

        write_delim(result,paste0("results_",FILE,"/haplotype_",filename,"_",FILE,"_",REGION,"_dominant.csv"), delim = ",", na = "NA", append = FALSE, col_names = TRUE, escape = "double")

}

haplotype_glm("DQA1|DQB1")
haplotype_glm("DQA1|DRB1")
haplotype_glm("DQB1|DRB1")
haplotype_glm("DQA1|DQB1|DRB1")

