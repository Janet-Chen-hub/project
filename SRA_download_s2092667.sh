#!/bin/sh 
#load the modules

module load igmm/apps/sratoolkit/2.8.2-1
module load igmm/apps/FastQC/0.11.4
module load roslin/subread/1.5.0-p3


### sra--------->fastq
#   everything is by default
#   for dbgap data all the work should be done in the working space like ~/ncbi/dbGAP_prj####
vdb-config --import prj_20892.ngc


mkdir QC
#subread-buildindex -o hg_index /exports/eddie/scratch/s2092667/Homo_sapiens.GRCh38.dna.primary_assembly.fa


while read acc;do
		echo "downloading sra files"
        prefetch $acc
        echo "converting ${acc} file to fastq"
        #paired-end
        cd /exports/eddie/scratch/s2092667/ncbi/dbGaP_20892/dbGaP-20892/sra
        fastq-dump --split-3 -gzip $acc.sra
        rm $acc.sra
        
        file1=$(ls -1|grep '_1.fastq.gz'|grep $acc)
		file2=$(ls -1|grep '_2.fastq.gz'|grep $acc)

        echo "fastqc started at $(date)"
        fastqc -o QCLIVER $file1
        fastqc -o QCLIVER $file2
        

done < /exports/eddie/scratch/s2092667/acc.txt


### quality control

mkdir QC
echo "fastqc started at $(date)"
fastqc -o QC *gz
echo "fastqc finished at $(date)"
## multiqc installation
##for convenience we could use multiqc to check the quality report 
## which is committed in local directory
conda create --name py3.7 python=3.7
conda activate py3.7
conda install -c bioconda -c conda-forge multiqc
cd QC
multiqc .
cd ..

#if it's on the eddie we need to first set the path and install it but it doesn't work now the roufh code is like
VERSION=0.7
firstPATH=/exports/eddie/scratch/s2092667/MultiQC
mkdir $firstPATH
INST=$firstPATH/$VERSION
#statements  
module load python/3.4.3
mkdir $INST

export PYTHONPATH=$INST/lib/python2.7/site-packages 
pip install --install-option="--prefix=$INST" multiqc

set components [ file split [ module-info name ] ]
set version [ lindex $components 1 ]
set modroot /path/to/software/multiqc/$version
proc ModulesHelp { } {
    global version modroot
    puts stderr "\tMultiQC - use MultiQC $version"
    puts stderr "\n\tVersion $version\n"
}
module-whatis   "Loads MultiQC environment."
# load required modules
module load python/2.7.6
# only one version at a time
conflict multiqc
# Make the directories available
prepend-path    PATH        $modroot/bin
prepend-path    PYTHONPATH  $modroot/lib/python2.7/site-packages



### subread
echo "subread started at $(date)"
##  build index
subread-buildindex -o hg_index Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

##  alignment and the counts 
mkdir csvlist
while read acc; do
	file1=$(ls -1|grep '_1.fastq.gz'|grep $acc)
	file2=$(ls -1|grep '_2.fastq.gz'|grep $acc)
	
	subread-align -t 1 -T 5 -i hg_index -r $file1 -R $file2 -o $acc.bam
	featureCounts -p -T 5 -t exon -g gene_id -a Homo_sapiens.GRCh38.98.gtf -o $acc.txt $acc.bam
    #remove the first two lines in the text files
    sed -i '1d' $acc.txt
    sed -i '1d' $acc.txt
    #count the mapped reads
    sum_count=$(cat $acc.txt|awk '{sum+=$7};END {print sum}')
    #count the fpkm values according to the formula and output to the csv files
    if [[ "$sum_count" != 0 ]]
    then
    	cat $acc.txt|awk '{OFS="\t";
    	{ print $6,$7,($7*(10^9)/($6*"'$sum_count'"))}
    	}'> $acc.csv
    	sed -i $'1 i\\\"'$acc"_length"'"\t"'$acc"_count"'"\t"'$acc"_fpkm"'"' $acc.csv #add a new row (col_names)
    fi
  
done < acc.txt

cp ./*csv ./csvlist
featureCounts -p -T 5 -t exon -g gene_id -a Homo_sapiens.GRCh38.98.gtf -o counts.txt ./*bam
echo "subread finished at $(date)"

