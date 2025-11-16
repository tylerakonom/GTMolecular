indirectory=/scratch/alpine/tyak9569/GTM/bqsr/

for pathandfilename in `ls ${indirectory}*_recal.bam`; do
name=`basename $pathandfilename _recal.bam`
echo $pathandfilename
echo $name
sbatch --export=filename=$name /scratch/alpine/tyak9569/GTM/bqsr/caller.sh 
done