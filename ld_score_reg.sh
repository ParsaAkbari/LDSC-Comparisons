#!/bin/bash
source ${OUT_DIR}/python_env/miniconda/bin/activate cond_analysis
this_pheno=$1
readarray -t PHENO_NAMES_ARRAY < ${PHENO_NAMES}

rg_list=$OUT_DIR/ldscore_regression/input/${this_pheno}.sumstats.gz
## now delete first element of the aray
## indexing of array when deleting starts at 1, but when picking the value (see this_pheno) it starts at 0
#unset PHENO_NAMES_ARRAY[${LSB_JOBINDEX}-1]

for pheno in ${PHENO_NAMES_ARRAY[@]}; do
    if [ "$pheno" == "$this_pheno" ]
    then
	continue
    fi
    rg_list=$rg_list,$OUT_DIR/ldscore_regression/input/${pheno}.sumstats.gz
done
#    --ref-ld-chr $SUPPORT_FILES/ldsc_files/eur_w_ld_chr/ \
#    --w-ld-chr $SUPPORT_FILES/ldsc_files/eur_w_ld_chr/ \

$SUPPORT_FILES/ldsc/ldsc.py \
    --rg ${rg_list} \
    --ref-ld-chr /lustre/scratch115/realdata/mdt2/projects/ukbiobank_t151/ldsc_herit/ldscores_precomp/eur_w_ld_chr/ \
    --w-ld-chr /lustre/scratch115/realdata/mdt2/projects/ukbiobank_t151/ldsc_herit/ldscores_precomp/eur_w_ld_chr/ \
    --out $OUT_DIR/ldscore_regression/output/ldreg_${this_pheno}

python_return_code=$?
if [ $python_return_code -ne 0 ]; then
    echo "Python script error"
    exit
fi

# we need to extract the portion of the log file which contains the correlation matrix
csplit --prefix ${this_pheno}xx $OUT_DIR/ldscore_regression/output/ldreg_${this_pheno}.log '/Summary of Genetic Correlation Results/' '/Analysis finished/'
rm ${this_pheno}xx00 ${this_pheno}xx02
mv ${this_pheno}xx01 $OUT_DIR/ldscore_regression/output/ldreg_${this_pheno}.txt

# delete the first line of the output file which is just a title that will confuse R
sed '1d' $OUT_DIR/ldscore_regression/output/ldreg_${this_pheno}.txt > $OUT_DIR/ldscore_regression/output/tempfile_${this_pheno}
mv $OUT_DIR/ldscore_regression/output/tempfile_${this_pheno} $OUT_DIR/ldscore_regression/output/ldreg_${this_pheno}.txt

# the phenotype names are the path to the input files each input file is one phenotype
# clean up the phenotype names to remove the part of the string which is the filepath
#sed -i 's@'$OUT_DIR'/ldscore_regression/input/@@g' $OUT_DIR/ldscore_regression/output/ldreg_${this_pheno}
sed -i 's@/lustre[^ ]*/input/@@g' $OUT_DIR/ldscore_regression/output/ldreg_${this_pheno}.txt
sed -i 's/.sumstats.gz//g' $OUT_DIR/ldscore_regression/output/ldreg_${this_pheno}.txt
