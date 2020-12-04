#perl and hmmer are needed to run itsx
module load perl
module load hmmer/3.1-snap20121016

#download itsx and extract it only if it ITSx_1.1.2 is not already installed in folder
FILE=ITSx_1.1.2
if [ ! -d "$FILE" ]; then
	echo "downloading ITSx"
	wget https://microbiology.se/sw/ITSx_1.1.2.tar.gz
	tar -xvf ITSx_1.1.2.tar.gz
else
	echo "ITSx exists, skipping download."
fi

#move the output of count.seqs to the itsx folder and then enter the folder
echo "Copying mothur output1 to ITSx folder"
cp $INTERPREFIX.fasta ITSx_1.1.2/$INTERPREFIX.fasta
cd ITSx_1.1.2

echo "Running ITSx on $INTERPREFIX.fasta"
#this line runs ITSx and extracts ITS1 and ITS2, then saves them as new fasta files. It also makes a summary.
./ITSx -i $INTERPREFIX.fasta -o $ITSXOUT -t f --cpu $PROCESSORS --summary T --graphical F --fasta T --preserve T --save_regions ITS1,ITS2 --anchor HMM
