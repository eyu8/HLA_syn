# Load required software modules
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# Activate the HLA-TAPAS virtual environment
source ~/runs/eyu8/env/HLA-TAPAS/bin/activate

# Set the R library path for HLA-TAPAS
export R_LIBS_SITE="/home/eyu8/runs/eyu8/library/HLA-TAPAS"

# Run HLAassoc with the logistic regression model
python -m HLAassoc LOGISTIC \
    --vcf ~/runs/eyu8/data/HLA_typing/TOPMED_HLA/iRBD.vcf \
    --out ~/runs/eyu8/data/HLA_typing/TOPMED_additive/TOPMED_additive \
    --pheno ~/runs/eyu8/data/HLA_typing/TOPMED_additive/iRBD_pheno.txt \
    --pheno-name pheno \
    --covar ~/runs/eyu8/data/HLA_typing/TOPMED_additive/iRBD_covar.txt

# Run HLAassoc with the omnibus logistic regression model using Slurm to allocate resources
srun -c 40 --mem=160g -t 1:0:0 python -m HLAassoc OMNIBUS_LOGISTIC \
    --vcf /home/eyu8/runs/eyu8/data/HLA_typing/TOPMED_additive/iRBD.vcf \
    --bim resources/wgsMHC.4digit.bglv4.bim \
    --fam ~/runs/eyu8/data/HLA_typing/TOPMED_additive/iRBD.fam \
    --out ~/runs/eyu8/data/HLA_typing/TOPMED_additive/TOPMED_additive_OMNIBUS \
    --aa-only \
    --pheno ~/runs/eyu8/data/HLA_typing/TOPMED_additive/iRBD_OMNIBUS_pheno.txt \
    --covar ~/runs/eyu8/data/HLA_typing/TOPMED_additive/iRBD_covar.txt \
    --maf-threshold 0
