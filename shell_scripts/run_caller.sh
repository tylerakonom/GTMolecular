indirectory=/scratch/alpine/tyak9569/GTM/

for pathandfilename in `ls ${indirectory}*.fastq`; do
name=`basename $pathandfilename .fastq`
echo $pathandfilename
echo $name
sbatch --export=filename=$name /scratch/alpine/tyak9569/GTM/alignReads.sh 
done