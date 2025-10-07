#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20gb		
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch
set -eou pipefail

# ENVIRONMENT SETUP
module purge
module load anaconda/25.3
module load seqkit/2.6
conda activate $NRDSTOR/repeat_annnotation

# Need conda environments for the various programs 
# reapeat_annotation == repeat masker & modeler 
# terl_env == TERL
# deepTE_env == DeepTE 

# Otherwise Copy my YAML file 

# [Repos] 
# DeepTE (https://github.com/LiLabAtVT/DeepTE?tab=readme-ov-file)

# Terl (https://github.com/muriloHoracio/TERL/tree/master)

# Repeat Annotation Pipeline (https://github.com/pedronachtigall/Repeat-annotation-pipeline)

# Michael Grapin @ Moore Lab UNL 
# ====================================================================================================== #
# Purpose: Annotate Repeat Regions of a Genome Assembly
# ====================================================================================================== #
WKDIR=

FASTA=

SPECIES=""

# DeepTE 
DEEPTE=/GitRepos/DeepTE/DeepTE.py
MODEL_TE=/home/moorelab/mgrapin2/SCRIPTS/HiFi_Pipeline/ANNOTATION/repeatMM_scripts/Metazoans_model/
# wget https://raw.githubusercontent.com/pedronachtigall/Repeat-annotation-pipeline/refs/heads/main/scripts/CleanDeepTEheader.py
CLEAN=/repeatMM_scripts/CleanDeepTEheader.py

# Terl 
TERL=/GitRepos/TERL/terl_test.py
MODEL_TERL=/repeatMM_scripts/DS3
# wget https://raw.githubusercontent.com/pedronachtigall/Repeat-annotation-pipeline/refs/heads/main/scripts/FilterTERL.py
FILTER=/repeatMM_scripts/FilterTERL.py


# Optional 
#wget https://raw.githubusercontent.com/pedronachtigall/Repeat-annotation-pipeline/refs/heads/main/scripts/AdjustingGFF_RM.py
GFF=/repeatMM_scripts/AdjustingGFF_RM.py

###########################################################################################################
##                                           Commamds                                                    ##
###########################################################################################################
mkdir -p $WKDIR/RepeatMM $WKDIR/std.err $WKDIR/std.out $WKDIR/RepeatMM/ML $WKDIR/RepeatMM/RepeatMasker
cd $WKDIR/RepeatMM

BuildDatabase -name $SPECIES $FASTA

RepeatModeler -threads $SLURM_CPUS_PER_TASK -database $SPECIES -numAddlRounds 2

echo "Done Running Repeat Modeler"


# Make a Variable for the Repeats.fa 
REPEAT=$WKDIR/RepeatMM/${SPECIES}-families.fa
cd $WKDIR/RepeatMM/ML

# add a prefix
cat $REPEAT | seqkit fx2tab | awk -v sp="$SPECIES" '{ print sp"_"$0 }' | seqkit tab2fx > "${SPECIES}_families.prefix.fa"

# separate files into known and unknown
cat ${SPECIES}_families.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > ${SPECIES}_families.prefix.fa.known
cat ${SPECIES}_families.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > ${SPECIES}_families.prefix.fa.unknown

# known set is ready to go
awk '/^>/{print $1; next}{print}' ${SPECIES}_families.prefix.fa.known > ${SPECIES}_families.prefix.fa.known.FINAL

echo "Running DeepTE..."
conda deactivate 
conda activate $NRDSTOR/deepTE_env
# Classifiy Unknown
python $DEEPTE -d deepTE_out -o deepTE_out -i ${SPECIES}_families.prefix.fa.unknown -sp M -m_dir $MODEL_TE
python $CLEAN deepTE_out/opt_DeepTE.fasta ${SPECIES}_families.prefix.fa.unknown.DeepTE

# Reporting
echo "Orginal Unknown: $(cat ${SPECIES}_families.prefix.fa.unknown | seqkit fx2tab | grep "Unknown" | wc -l)"
echo "Unknown After DeepTE: $(cat ${SPECIES}_families.prefix.fa.unknown.DeepTE | seqkit fx2tab | grep "Unknown" | wc -l)"

echo "Running TERL..."
conda deactivate 
conda activate $NRDSTOR/terl_env
python $TERL -m $MODEL_TERL -f ${SPECIES}_families.prefix.fa.unknown.DeepTE
mv TERL* ${SPECIES}_families.prefix.fa.unknown.DeepTE.TERL

conda deactivate 
conda activate $NRDSTOR/prepocess_repeats
echo "TERL Classification:"
python $FILTER ${SPECIES}_families.prefix.fa.unknown.DeepTE.TERL ${SPECIES}_families.prefix.fa.unknown.DeepTE ${SPECIES}_families.prefix.fa.unknown.FINAL

# Reporting 
echo "Unknown After TERL: $(cat ${SPECIES}_families.prefix.fa.unknown.FINAL | seqkit fx2tab | grep "#Unknown" | wc -l)"

# Run RepeatMasker Again to generate full annotation

# concatenate known and unknown
cat ${SPECIES}_families.prefix.fa.known.FINAL ${SPECIES}_families.prefix.fa.unknown.FINAL > ${SPECIES}_families.prefix.fa.known.unknown.FINAL

echo "Final De Novo Custom Repeat Library: $(grep -c ">" ${SPECIES}_families.prefix.fa.known.unknown.FINAL) Sequences" 


conda deactivate 
conda activate $NRDSTOR/repeat_annnotation
cd $WKDIR/RepeatMM/RepeatMasker
#run repeatmasker
RepeatMasker -pa $SLURM_CPUS_PER_TASK -gff -s -a -inv -no_is -norna -xsmall -nolow -div 40 -lib $WKDIR/RepeatMM/ML/${SPECIES}_families.prefix.fa.known.unknown.FINAL -cutoff 225 $FASTA

exit 

#OPTIONAL - calculate the kimura distance
calcDivergenceFromAlign.pl -s $SPECIES.fa.align.divsum $SPECIES.fa.align
tail -n 72 $SPECIES.fa.align.divsum > $SPECIES.fa.Kimura.distance

#OPTIONAL - generate a modified GFF file, which I like more than the default generated by RepeatMasker.
python $GFF $SPECIES.fa.out.gff ${SPECIES}_families.prefix.fa.known.unknown.FINAL $SPECIES.fa.out.adjusted.gff



echo "######################################################################"
echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"

echo "Job finished at: $(date)"
