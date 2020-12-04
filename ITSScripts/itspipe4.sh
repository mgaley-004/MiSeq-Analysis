#now calculate percent identity of the sequences using python
#to use this script, supply in order the output name you would like, the reference fasta file, the taxonomy file that mothur created for sequences (not the cons.taxonomy file for otus and not the equalized taxonomy file), and the last fasta file that mothur created. This script should fill those in for you if you've changed the variables at the top of the script.
PYSCRIPT2=computePercentIdentityParallel.py
if [ ! -f "$PYSCRIPT2" ]; then
	echo "downloading computePercentIdentityParallel.py"
	wget https://github.com/mgaley-004/MiSeq-Analysis/blob/main/ITSScripts/computePercentIdentityParallel.py
else
	echo "computePercentIdentityParallel.py exists, skipping download."
fi

echo "Computing percent identity"

python computePercentIdentityParallel.py $FINALFILENAME $REFERENCEFASTA $FINALMOTHUR.taxonomy $FINALMOTHUR.fasta
#from here we will use some sort of thresholding to assign confidence in taxonomic IDs.