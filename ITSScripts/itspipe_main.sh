#~/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=galey004@umn.edu
#SBATCH -p small
#SBATCH -o %j.out
#SBATCH -e %j.err

TERM=linux
export TERM
export FASTQLOC="/panfs/roc/groups/2/chun0157/data_release/umgc/miseq/200413_M00262_0047_000000000-J2LDF/Chun_Project_005_ITS_Redemult"
export SITES="*NP*.gz *BRL*.gz"
export PREFIX=its19
export INTERPREFIX=mothur1
export ITSXOUT=itsxo
export PROCESSORS=$(nproc)
export WDIR="/scratch.global/galey_ex"
export REFERENCEFASTA=/panfs/roc/groups/2/chun0157/shared/UNITE.ITS1.good.fasta
export REFERENCETAX=/panfs/roc/groups/2/chun0157/shared/ITS1.UNITE.good.full.taxonomy
export FINALMOTHUR=mothur2
export FINALFILENAME=pipeline_output.csv

if [ ! -d "$WDIR" ]; then
	echo "creating directory $WDIR"
	mkdir $WDIR
else
	echo "$WDIR exists. Make sure this directory is empty before beginning."
fi

cd $WDIR
BATCH1=its_part1.batch
BATCH2=its_part2.batch
if [ ! -f "$BATCH1" ]; then
	echo "downloading $BATCH1"
	wget https://raw.githubusercontent.com/mgaley-004/MiSeq-Analysis/main/ITSScripts/its_part1.batch
fi
if [ ! -f "$BATCH2" ]; then
	echo "downloading $BATCH2"
	wget https://raw.githubusercontent.com/mgaley-004/MiSeq-Analysis/main/ITSScripts/its_part2.batch
fi

#install the latest copy of mothur
wget https://github.com/mothur/mothur/releases/download/v1.44.3/Mothur.linux.zip
unzip Mothur.linux.zip

#If this is not the first time you have run this script, change the lines below by removing dependencies for steps already completed and then comment out the step lines.

job1=$(sbatch --parsable ~/itspipe1.sh)
job2=$(sbatch --parsable --dependency=afterok:${job1} ~/itspipe2.sh)
job3=$(sbatch --parsable --dependency=afterok:${job2} ~/itspipe3.sh)
job4=$(sbatch --parsable --dependency=afterok:${job3} ~/itspipe4.sh)

#Step 1 moves the files defined in SITES to WDIR, moves a copy of mothur to WDIR, and runs the first mothur batch script.
srun --nodes 1 --ntasks 24 job1

#Step 2 runs ITSx to identify ITS1 sequecnes from mothur output. Will not run if step 1 did not complete succesfully.
srun --nodes 1 --ntasks 24 job2

#Step 3 runs the second mothur batch script using output of step 2.
srun --nodes 1 --ntasks 24 job2

#Step 4 computes percent identity.
srun --nodes 1 --ntasks 24 job2