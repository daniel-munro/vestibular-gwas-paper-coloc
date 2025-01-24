# Extract list of tissues to run fastENLOC for each one
zcat data/eqtl/gtex_v8.eqtl_annot_rsid.vcf.gz \
    | cut -f6 \
    | sed 's/|/\n/g' \
    | cut -d'@' -f2 \
    | cut -d'=' -f1 \
    | sort \
    | uniq \
    > data/eqtl/tissues.txt

#####################
## Vestibular GWAS ##
#####################

# Prepare GWAS PIP (posterior inclusion probabilities) from summary stats
# Get z-scores from slope + slope SE:
python3 scripts/torus_input.py \
    --gwas data/gwas/EUR_balance180.gz \
    --format mvp \
    --ld data/reference/fourier_ls-all.EUR.hg19.bed \
    --out data/gwas/EUR_balance180.zval.gz
# Convert z-scores to PIPs:
# Segmentation fault probably means a file name is invalid
~/tools/torus/src/torus \
    -d data/gwas/EUR_balance180.zval.gz \
    --load_zval \
    -dump_pip data/gwas/EUR_balance180.gwas.pip
# Without -f I get "gzip: data/fastenloc/EUR_balance180.gwas.pip has 1 other link  -- unchanged"
gzip -f data/gwas/EUR_balance180.gwas.pip

# Run fastENLOC with precomputed GTEx inputs
# Using rsID version since I believe GWAS used GRCh37 and eQTL annotations are GTCh38
parallel -j4 \
    < data/eqtl/tissues.txt \
    ~/tools/fastenloc/src/fastenloc \
    -eqtl data/eqtl/gtex_v8.eqtl_annot_rsid.vcf.gz \
    -gwas data/gwas/EUR_balance180.gwas.pip.gz \
    -tissue {} \
    -thread 2 \
    -prefix data/fastenloc/EUR_balance180/EUR_balance180.{}

##########################
## GTEx Harmonized GWAS ##
##########################

# Prepare GWAS PIP (posterior inclusion probabilities) from summary stats
# Get z-scores from slope + slope SE:
python3 scripts/torus_input.py \
    --gwas data/gwas_gtex/imputed_UKB_50_Standing_height.txt.gz \
    --format gtex \
    --ld data/reference/eur_ld.hg38.bed \
    --out data/gwas_gtex/imputed_UKB_50_Standing_height.zval.gz
# Convert z-scores to PIPs:
# Segmentation fault probably means a file name is invalid
~/tools/torus/src/torus \
    -d data/gwas_gtex/imputed_UKB_50_Standing_height.zval.gz \
    --load_zval \
    -dump_pip data/gwas_gtex/imputed_UKB_50_Standing_height.gwas.pip
# Without -f I get "gzip: data/fastenloc/EUR_balance180.gwas.pip has 1 other link  -- unchanged"
gzip -f data/gwas_gtex/imputed_UKB_50_Standing_height.gwas.pip

# Run fastENLOC with precomputed GTEx inputs
# Using rsID version since I believe GWAS used GRCh37 and eQTL annotations are GTCh38
mkdir -p data/fastenloc/imputed_UKB_50_Standing_height
parallel -j4 \
    < data/eqtl/tissues.txt \
    ~/tools/fastenloc/src/fastenloc \
    -eqtl data/eqtl/gtex_v8.eqtl_annot_rsid.vcf.gz \
    -gwas data/gwas_gtex/imputed_UKB_50_Standing_height.gwas.pip.gz \
    -tissue {} \
    -thread 2 \
    -prefix data/fastenloc/imputed_UKB_50_Standing_height/imputed_UKB_50_Standing_height.{}
