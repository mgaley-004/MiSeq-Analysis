#to resume mothur, we have to remove sequences that didn't contain ITS1 from the count file. We will create a .accnos file using the original input and the ITS1 fasta files

#load python
module load python
#run the create_accnos python script
PYSCRIPT1=create_accnos.py
if [ ! -f "$PYSCRIPT1" ]; then
	echo "downloading create_accnos.py"
	wget https://github.com/mgaley-004/MiSeq-Analysis/blob/main/ITSScripts/create_accnos.py
else
	echo "create_accnos.py exists, skipping download."
fi

echo "Running create_accnos.py"
python create_accnos.py $ITSXOUT.ITS1.fasta $INTERPREFIX.fasta nontargets.accnos
cp nontargets.accnos $WDIR/nontargets.accnos
cp $ITSXOUT.ITS1.fasta $WDIR/$ITSXOUT.ITS1.fasta

echo "Resuming mothur to classify sequences."
#resume mothur
mothur/./mothur "its_part2.batch"