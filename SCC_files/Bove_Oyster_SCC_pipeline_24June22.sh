## Eastern Oyster OA and boring sponge project
# Pipeline by Colleen B Bove
#------------------------------

### Move files from TUFC to SCC

## Log on to your scc account using ssh: ssh bovec@scc1.bu.edu
# log on to TUCF’s ftp 
# adding -i will get rid of prompts later on
ftp -i 

## you will be prompted for username and password
# Username: 
# Password: 
## you might have to type the following to get ftp to prompt you fo username
user

# NOTE: tab fills are deactivated in ftp

## to pull multiple files to remote server use
mget *fastq bovec@scc1.bu.edu:/projectnb/davieslab/bove/raw_seqs/oysters_raw

# if you get an error like WARNING! 1433673 bare linefeeds received in ASCII mode 
# this means there was a problem entering binary mode and even if you successfully transfer files, it will screw up any zipped files. 
# to get around this before doing mget do 
binary 
# then do mget as described above

# if you get an error like 426 Failure writing network stream
# this is likely because the remote server took too long to download the file, you can try #when you have better internet, or try downloading one file at a time. 

## -----
# Copy the fastq files to working folder
cp /projectnb/davieslab/bove/raw_seqs/oysters_raw/*.gz /projectnb/davieslab/bove/oysters/fastq_files

# unzip all files (I didn't submit as a job here)
gunzip *.gz

# just for sanity check, make sure a file looks correct:
head -50 19650-Bove_oyster_S22_R1_001.fastq  # looks good!
# note that every read has four lines, the ID line starts with @D00780



## Be sure to do TUFT’s recommended integrity check after unzipping your files!! 
# compare md5sum outputs by first generating one
md5sum *.gz > checklist.chk

# and then check that they match
diff checklist.chk report_files/md5sum.txt

# check the number of files
ls | grep .fastq | wc -l


## remove unneeded text from file names (do this for .gz files)
for f in *.fastq; do
mv -- "$f" "${f/-Bove_oyster/}"
done

for f in *.fastq; do
mv -- "$f" "${f/_R1_001/}"
done




### Adaptor and quality trimming:

# load fastx_toolkit
module load fastx-toolkit/0.0.14
module load perl/5.28.1

# creating and launching the cleaning process for all files in the same time (job called 'clean')

for file in *.fastq
do echo "perl tagseq_clipper.pl ${file} | fastx_clipper -a AAAAAAAA -l 20 | fastx_clipper -a AGATCGGAAG -l 20 | fastq_quality_filter -q 20 -p 90 -Q33 >${file/.fastq/}.trim" >>clean
done


#-------------
#	$ fastx_clipper -h
#	usage: fastx_clipper [-h] [-a ADAPTER] [-D] [-l N] [-n] [-d N] [-c] [-C] [-o] [-v] [-z] [-i INFILE] [-o OUTFILE]
#
#	version 0.0.6
#	   [-a ADAPTER] = ADAPTER string. default is CCTTAAGG (dummy adapter).
#	   [-l N]       = discard sequences shorter than N nucleotides. default is 5.
#	   [-d N]       = Keep the adapter and N bases after it.
#	   [-c]         = Discard non-clipped sequences (i.e. - keep only sequences which contained the adapter).
#	   [-C]         = Discard clipped sequences (i.e. - keep only sequences which did not contained the adapter).
#	   [-k]         = Report Adapter-Only sequences.
#			  (using '-d 0' is the same as not using '-d' at all. which is the default).
#	   [-n]         = keep sequences with unknown (N) nucleotides. default is to discard such sequences.
#	   [-v]         = Verbose - report number of sequences.
#			  If [-o] is specified,  report will be printed to STDOUT.
#			  If [-o] is not specified (and output goes to STDOUT),
#			  report will be printed to STDERR.
#	   [-z]         = Compress output with GZIP.
#	   [-D]		= DEBUG output.
#	   [-i INFILE]  = FASTA/Q input file. default is STDIN.
#	   [-h]         = This helpful help screen.
#	   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.

#------

#$ fastq_quality_filter -h
#	usage: fastq_quality_filter [-h] [-v] [-q N] [-p N] [-z] [-i INFILE] [-o OUTFILE]
#
#	version 0.0.6
#	   [-h]         = This helpful help screen.
#	   [-q N]       = Minimum quality score to keep.
#	   [-p N]       = Minimum percent of bases that must have [-q] quality.
#	   [-z]         = Compress output with GZIP.
#	   [-i INFILE]  = FASTA/Q input file. default is STDIN.
#	   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.
#	   [-v]         = Verbose - report number of sequences.
#			  If [-o] is specified,  report will be printed to STDOUT.
#			  If [-o] is not specified (and output goes to STDOUT),
#			  report will be printed to STDERR.

#-------------


# add the bash heading and run as shell script
# ---
#!/bin/bash -l

# Give job a name
#$ -N clean_hyp

# Combine output and error files into a single file
#$ -j y
# ---


# submit job
qsub clean.sh

# check the status of my jobs running
qstat -u bovec


# It is complete! I got a bunch of .trim files that are non-empty! 
# but did the trimming really work? 
# Use the same one-liner as before on the trimmed file to see if it is different
# from the raw one that you looked at before:
head -50 19637_S6.trim | grep -E '^[ACGT]+$'


# ------------
### What the cleaning/trimming does:

## tagseq_clipper.pl:
# "Clips 5'-leader off Illumina fastq reads in RNA-seq
# Removes duplicated reads sharing the same degenerate header and
# the first 20 bases of the sequence (reads containing N bases in this
# region are discarded, too)"

## fastx_clipper:
# [-a ADAPTER] = ADAPTER string. default is CCTTAAGG (dummy adapter). ----> I use AAAAAAAA and AGATCGGAAG with
# [-l N]       = discard sequences shorter than N nucleotides. default is 5. ----> I used 20 for both adapters

## fastq_quality_filter:    fastq_quality_filter -q 20 -p 90 -Q33
# [-q N]       = Minimum quality score to keep. ----> I used 20 nucleotides
# [-p N]       = Minimum percent of bases that must have [-q] quality. ----> I used 90 % here
# ------------





### Mapping to C. virginica genome:
## genome was accessed from NCBI: https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica
# Get the C. virginica genome direct from NCBI:
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/* /projectnb/davieslab/bove/oysters/genome_files 



## Need to get gffread from source since it is not on the SCC
# Download from: http://ccb.jhu.edu/software/stringtie/gff.shtml
# download the source package version (this was successful) (location: /projectnb/davieslab/bove/gffread-0.12.7)
# followed the git installation instructions (https://github.com/gpertea/gffread)
cd gffread-0.12.7/
make release

# by running the following code, we are now able to run featureCounts from outside this directory
export PATH=/projectnb/davieslab/bove/gffread-0.12.7/gffread-0.12.7:$PATH




# Make shell script to convert from .gff to .gtf
nano gff2gtf.sh

# ---
#!/bin/bash -l

#$ -j y

gffread GCF_002022765.2_C_virginica-3.0_genomic.gff -T -o Cvir_genome.gtf
# ---

# We now have a non-empty genome file in gtf format: Cvir_genome.gtf



## Make a STAR index with CVIR annotation files

# Make a shell script to run:

nano STAR_genomeCreate.sh



# ---
#!#/bin/bash
#$ -j y

echo "Please put in the base directory:"
read base
echo "Please put in the output folder name:"
read output

echo "Outputs saving to : " $base$output

if [ -d "$base$output" ]; then
	echo "Directory Already Exists, please rerun with unique output directory"
	exit 1
else
	echo "Directory Created"
	mkdir "$base$output"
fi

echo "Select genome file (.fna format, should include entire path)"
read genome

echo "Select gene annotation file (.gtf, should includ entire path)"
read gene_annotation


STAR --runThreadN 32 \
--runMode genomeGenerate \
--genomeDir $base$output \
--genomeFastaFiles $genome \
--sjdbGTFfile $gene_annotation \
--genomeSAindexNbases 13

# ---

# Initially ran it without the '--genomeSAindexNbases 13' line and got the following warming (output in star_ref):
# !!!!! WARNING: --genomeSAindexNbases 14 is too large for the genome size=684741128, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 13
# ... so reran the script and saved it to: STAR_ref 

## Originally mapped to the NCBI genome (GCF_002022765.2_C_virginica-3.0_genomic.fna), however, this version of the genome contained a lot of haplotigs.
# Because of this, using a new masked genome from the Puritz lab.
# GitHub repo: https://github.com/The-Eastern-Oyster-Genome-Project/2022_Eastern_Oyster_Haplotig_Masked_Genome/tree/main/Haplotig_Masking/Output/Masked_Genome
# biorxiv manuscript:  https://www.biorxiv.org/content/10.1101/2022.08.29.505626v1

# move the masked genome to the SCC
scp -r CVIR_reference_masked.fasta bovec@scc1.bu.edu:/projectnb/davieslab/bove/oysters/genome_files



# load the STAR module
module load star/2.7.9a

# run the STAR script
source STAR_genomeCreate.sh


### Input requests from running the above script
# Please put in the base directory:
# /projectnb/davieslab/bove/oysters/genome_files/
# 
# Please put in the output folder name:
# STAR_ref_mask
# 
# Outputs saving to :  /projectnb/davieslab/bove/oysters/genome_files/gffread-0.12.7/STAR_ref
# Directory Created
# 
# Select genome file (.fna format, should include entire path)
# /projectnb/davieslab/bove/oysters/genome_files/CVIR_reference_masked.fasta
# 
# Select gene annotation file (.gtf, should includ entire path)
# /projectnb/davieslab/bove/oysters/genome_files/GCF_002022765.2_C_virginica-3.0_genomic.gtf


## Output from the successful run:
# 	STAR --runThreadN 32 --runMode genomeGenerate --genomeDir /projectnb/davieslab/bove/oysters/genome_files/STAR_ref_mask --genomeFastaFiles /projectnb/davieslab/bove/oysters/genome_files/CVIR_reference_masked.fasta --sjdbGTFfile /projectnb/davieslab/bove/oysters/genome_files/GCF_002022765.2_C_virginica-3.0_genomic.gtf --genomeSAindexNbases 13
# 	STAR version: 2.7.9a   compiled: 2021-07-19T14:08:30-0400 scc-bc1:/share/pkg.7/star/2.7.9a/src/STAR-2.7.9a/source
# Sep 23 09:17:01 ..... started STAR run
# Sep 23 09:17:01 ... starting to generate Genome files
# Sep 23 09:17:07 ..... processing annotations GTF
# Sep 23 09:17:18 ... starting to sort Suffix Array. This may take a long time...
# Sep 23 09:17:20 ... sorting Suffix Array chunks and saving them to disk...
# Sep 23 09:18:40 ... loading chunks from disk, packing SA...
# Sep 23 09:18:54 ... finished generating suffix array
# Sep 23 09:18:54 ... generating Suffix Array index
# Sep 23 09:19:49 ... completed Suffix Array index
# Sep 23 09:19:50 ..... inserting junctions into the genome indices
# Sep 23 09:21:03 ... writing Genome to disk ...
# Sep 23 09:21:04 ... writing Suffix Array to disk ...
# Sep 23 09:21:12 ... writing SAindex to disk
# Sep 23 09:21:13 ..... finished successfully



## Making a new directory within 'oysters' for the output
# using mapped_files as of 23 June 2022 because of issues with the previous run (old ones are in 'old_mapping')
mkdir mapped_files



## STAR Mapping Test
# Going to test the commands to make sure they work with the following sample: 19013_S27.fastq


nano STAR_mapping.sh

#----
#!/bin/bash -l

STAR --genomeDir /projectnb/davieslab/bove/oysters/genome_files/STAR_ref_mask \
      --readFilesIn 19013_S27.fastq \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outSAMunmapped Within \
      --outFileNamePrefix ../mapped_files/19013_S27      
#----



## This worked! Next step will be to try to run this to loop through all 
# quick modification to the script to be able to loop through samples and run all as job

nano STAR_mapping.sh

#----
#!/bin/bash -l

for fastq in /projectnb/davieslab/bove/oysters/trim_files/*trim
do echo $fastq
$RUN STAR --genomeDir /projectnb/davieslab/bove/oysters/genome_files/STAR_ref_mask \
      --readFilesIn $fastq \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outSAMunmapped Within \
      --outFileNamePrefix ../mapped_files/$(basename $fastq trim)
done
#----

# submit this job
qsub nano STAR_mapping.sh

# and see if it runs!
qstat -u bovec

# we should also see a number of files per sample name being created within the 'star_mapped' directory


# use this to view mapped read % for each sample
grep "Uniquely mapped reads %" *Log.final.out

# 19013_S27.Log.final.out:                        Uniquely mapped reads % |	78.81%
# 19041_S8.Log.final.out:                        Uniquely mapped reads % |	77.41%
# 19076_S18.Log.final.out:                        Uniquely mapped reads % |	77.26%
# 19086_S16.Log.final.out:                        Uniquely mapped reads % |	75.73%
# 19092_S11.Log.final.out:                        Uniquely mapped reads % |	76.86%
# 19103_S5.Log.final.out:                        Uniquely mapped reads % |	75.60%
# 19158_S29.Log.final.out:                        Uniquely mapped reads % |	77.29%
# 19197_S21.Log.final.out:                        Uniquely mapped reads % |	79.24%
# 19227_S4.Log.final.out:                        Uniquely mapped reads % |	76.85%
# 19254_S12.Log.final.out:                        Uniquely mapped reads % |	78.05%
# 19264_S19.Log.final.out:                        Uniquely mapped reads % |	75.52%
# 19274_S7.Log.final.out:                        Uniquely mapped reads % |	76.09%
# 19304_S31.Log.final.out:                        Uniquely mapped reads % |	77.34%
# 19355_S1.Log.final.out:                        Uniquely mapped reads % |	71.65%
# 19370_S30.Log.final.out:                        Uniquely mapped reads % |	75.80%
# 19372_S17.Log.final.out:                        Uniquely mapped reads % |	77.47%
# 19378_S23.Log.final.out:                        Uniquely mapped reads % |	77.62%
# 19382_S24.Log.final.out:                        Uniquely mapped reads % |	70.80%
# 19388_S20.Log.final.out:                        Uniquely mapped reads % |	76.98%
# 19441_S15.Log.final.out:                        Uniquely mapped reads % |	77.83%
# 19467_S14.Log.final.out:                        Uniquely mapped reads % |	76.87%
# 19472_S26.Log.final.out:                        Uniquely mapped reads % |	72.84%
# 19496_S9.Log.final.out:                        Uniquely mapped reads % |	76.64%
# 19505_S25.Log.final.out:                        Uniquely mapped reads % |	73.71%
# 19516_S13.Log.final.out:                        Uniquely mapped reads % |	78.52%
# 19536_S3.Log.final.out:                        Uniquely mapped reads % |	76.47%
# 19549_S28.Log.final.out:                        Uniquely mapped reads % |	72.82%
# 19611_S2.Log.final.out:                        Uniquely mapped reads % |	73.21%
# 19637_S6.Log.final.out:                        Uniquely mapped reads % |	79.09%
# 19642_S10.Log.final.out:                        Uniquely mapped reads % |	77.56%
# 19650_S22.Log.final.out:                        Uniquely mapped reads % |	77.42%
# Undetermined_S0.Log.final.out:                        Uniquely mapped reads % |	78.15%



## Going to check that a bam file looks good using samtools:
module load samtools
samtools view -h 19650_S22.Aligned.sortedByCoord.out.bam | less

# this looks good so moving forward!



#------------------------------
## Generating read-counts-per gene
# this is a helpful walkthrough: https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/master/lessons/2day_rnaseq_workflow.md


## Need to add featureCounts to SCC:
# Download tar.gz from here: http://subread.sourceforge.net/
# move to the SCC 
scp -r subread-2.0.3-source.tar.gz bovec@scc1.bu.edu:/projectnb/davieslab/bove/subread-2.0.3

tar -xzf subread-2.0.3-source.tar.gz
cd subread-2.0.3-source/src
make -f Makefile.Linux

# It works, now this is the path needed for featureCounts:
/projectnb/davieslab/bove/subread-2.0.3/subread-2.0.3-source/bin

# by running the following code, we are now able to run featureCounts from outside this directory
export PATH=/projectnb/davieslab/bove/subread-2.0.3/subread-2.0.3-source/bin:$PATH



## Move all .bam files from 'mapped_files' to 'bam_files'
mv mapped_files/*.bam bam_files
cd bam_files


## Now that we have featureCounts sorted, lets give it a try with a single sample
# need to remove the emptyp gene_id from the GCF_002022765.2_C_virginica-3.0_genomic.gtf file for the following to run : GCF_002022765.2_C_virginica-3.0_genomic_UPDATE.gtf
grep 'gene_id ""' GCF_002022765.2_C_virginica-3.0_genomic.gtf 
grep -v 'gene_id ""' GCF_002022765.2_C_virginica-3.0_genomic.gtf > GCF_002022765.2_C_virginica-3.0_genomic_UPDATE.gtf


# make a new 'count_files' folder 
mkdir count_files

featureCounts -T 6 -s 1 \
  -a /projectnb/davieslab/bove/oysters/genome_files/GCF_002022765.2_C_virginica-3.0_genomic_UPDATE.gtf \
  -o /projectnb/davieslab/bove/oysters/count_files/CVIR_featurecounts_masked.txt \
  /projectnb/davieslab/bove/oysters/bam_files/*bam



# Now check that it looks okay (the first few columns will look a bit funky, but we will remove them!)
less /projectnb/davieslab/bove/oysters/count_files/CVIR_featurecounts_masked.txt



## Done! Just tidying up the 
# shorten the column name to sample name only
sed 's/.Aligned.sortedByCoord.out.bam//g' CVIR_featurecounts_masked.txt>CVIR_featurecounts_masked2.txt
sed 's:/projectnb/davieslab/bove/oysters/bam_files/::g' CVIR_featurecounts_masked2.txt>CVIR_featurecounts_23Sep22.txt

# remove the first line (it is not necessary)
sed -i '1d' CVIR_featurecounts_23Sep22.txt 




#------------------------------
## All done!

# scp CVIR_featurecounts_22Jun22.txt file to your computer
scp bovec@scc1.bu.edu:/projectnb/davieslab/bove/oysters/count_files/CVIR_featurecounts_23Sep22.txt /Users/Colleen/Dropbox/BU/DaviesLab/Oysters/SCC_pipeline_oyster



####### For downstream, we will need a couple other files
## Here is the link to those files: https://github.com/tejashree1modak/Cvir-cnv/tree/master/annotation
# we are downloading the following:
ref_annot_prot.tab # gene ID to GO annotation ID
ref_annot.txt # this is the gene ID to the full gene name
annotations.txt # GO annotation ID to GO terms