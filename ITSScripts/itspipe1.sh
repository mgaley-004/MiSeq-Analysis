#copy the needed files from the shared directory and the 
cp -r $MOTHURLOC $WDIR
cd $FASTQLOC
cp $SITES $WDIR

#change to scratch directory
cd $WDIR

#this first mothur batch script processes the fastq files and exports raw sequences
mothur/./mothur "its_part1.batch"