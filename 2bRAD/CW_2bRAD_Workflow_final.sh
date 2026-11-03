#This document details the data acquisition and analysis of siderastrea siderea and branching porites species
#The goal of this analysis is to look for cryptic lineages of both species in Curaçao from two distinct environments (bay and reef)

#Generally following Misha Matz's pipeline, found here: https://github.com/z0on/2bRAD_denovo

#Download scripts from https://github.com/z0on/2bRAD_denovo
#[haich@scc1 TVE_2bRAD]$ pwd
#/projectnb/davies-hb/hannah/TVE_Panama/TVE_2bRAD
/proj/kdcastil/users/mayapow/2bRAD


git clone https://github.com/z0on/2bRAD_denovo.git

# Make everything executable
chmod +x *.pl
chmod +x *.py
chmod +x *.R

# Install modules - we will need bowtie2, samtools, and picard. They are pre-installed as modules on TACC; you will have to install them if you don't have these modules on your cluster

module load perl
module load bowtie2
module load samtools
module load picard


# Working directory for this project is:
pwd
#/projectnb/davies-hb/hannah/TVE_2bRAD
#/proj/kdcastil/users/mayapow/2bRAD
/proj/kdcastil/users/mayapow/2bRAD/denovo_por #por just my samples
/proj/kdcastil/users/mayapow/2bRAD/denovo_sid #sid just my samples
/proj/kdcastil/users/mayapow/2bRAD/denovo_nosym #sid all caribbean samples (hannah & belize)

###MAYA STARTING HERE!!!####
#Deleted everything above here since it was done for me thanks Ally!!
#Download scripts from https://github.com/z0on/2bRAD_denovo
#[haich@scc1 TVE_2bRAD]$ pwd
#/projectnb/davies-hb/hannah/TVE_Panama/TVE_2bRAD
[mayapow@longleaf-login1 2bRAD]$ pwd
/proj/kdcastil/users/mayapow/2bRAD

git clone https://github.com/z0on/2bRAD_denovo.git

# Make everything executable (go to 2bRAD_denovo then do this)
cd 2bRAD_denovo
chmod +x *.pl
chmod +x *.py
chmod +x *.R

# Install modules - we will need bowtie2, samtools, and picard. They are pre-installed as modules on TACC; you will have to install them if you don't have these modules on your cluster
cd ..
module load perl
module load bowtie2
module load samtools
module load picard


# Working directory for this project is:
pwd
#/projectnb/davies-hb/hannah/TVE_2bRAD
/proj/kdcastil/users/mayapow/2bRAD

#------------------------------DE NOVO RAD BUSINESS

#get hannah files
#PRJNA833675
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA973074&result=read_run&fields=run_accession,sample_alias,fastq_ftp,fastq_md5&format=tsv" > filereport.tsv

cut -f3 filereport.tsv | tail -n +2 | tr ';' '\n' | sed 's|^|ftp://|' > urls.txt
wget -c -i urls.txt

#get karina files
#PRJNA833675
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA833675&result=read_run&fields=run_accession,sample_alias,fastq_ftp,fastq_md5&format=tsv" > filereport.tsv

cut -f3 filereport.tsv | tail -n +2 | tr ';' '\n' | sed 's|^|ftp://|' > urls.txt
wget -c -i urls.txt


#move to new folder:
mv /proj/kdcastil/users/mayapow/2bRAD/*.fastq.gz /proj/kdcastil/users/mayapow/2bRAD/denovo_sid_combo
scp -r 

# mv all *nosymbio* files to a new direcotry

#mkdir denovo_nosym
#mv *nosymbio* denovo_nosym/
mv maya-fastqs denovo_nosym

#pwd
#/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym
pwd
#/proj/kdcastil/users/mayapow/2bRAD/denovo_nosym
/proj/kdcastil/users/mayapow/2bRAD/denovo_por

# 'uniquing' ('stacking') individual fastq reads:

#ls *.nosymbio.fastq | perl -pe 's/^(.+)$/uniquerOne.pl $1 >$1\.uni/' >unii
#need to specify full filepath or it won't find uniquerOne.pl
ls *.fastq | perl -pe 's/^(.+)$/\/proj\/kdcastil\/users\/mayapow\/2bRAD\/2bRAD_denovo\/uniquerOne.pl $1 >$1\.uni/' >unii

head -20 unii #check out your list of files to make sure they have the unii tag

bash unii
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N unii # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

#../../2bRAD_denovo/uniquerOne.pl 1-MullenDavies_S1_ACCA.nosymbio.fastq >1-MullenDavies_S1_ACCA.nosymbio.fastq.uni
#../../2bRAD_denovo/uniquerOne.pl 1-MullenDavies_S1_AGAC.nosymbio.fastq >1-MullenDavies_S1_AGAC.nosymbio.fastq.uni
#../../2bRAD_denovo/uniquerOne.pl 1-MullenDavies_S1_AGTG.nosymbio.fastq >1-MullenDavies_S1_AGTG.nosymbio.fastq.uni
#../../2bRAD_denovo/uniquerOne.pl 1-MullenDavies_S1_CATC.nosymbio.fastq >1-MullenDavies_S1_CATC.nosymbio.fastq.uni

# uniquerOne.pl :
# Makes uniqued 2bRAD read file for a single fastq file. 
#This is analogous to making 'stacks' in STACKS. 
#The script records unique tag sequences, the total number of their appearances, and how many of those were in reverse-complement orientation. 
#(STACKS would consider reverse-complements as separate loci)

# merging uniqued files (set minInd to >10, or >10% of total number of samples, whatever is greater)
# since I have 61 samples I am keeping the minInd at 10

#I have 61 (Hannah)+ 24 (me) = 85
#for denovo_nosym full dataset I have 61 (Hannah) + 24 (me) +  (belize) = 85

#../../2bRAD_denovo/mergeUniq.pl uni minInd=10>all.uniq
/proj/kdcastil/users/mayapow/2bRAD/2bRAD_denovo/mergeUniq.pl uni minInd=10>all.uniq

# 313722 tags before tufts custom code
# 178617 tags after

#considering reads seen at least 5 times:
#47777868 reads processed maya
#39397105 reads processed maya porites
#8209907 reads processed maya siderastrea
#considering reads seen at least 5 times:
#30485376 reads processed maya sid combo hannah
#considering reads seen at least 5 times:
#70673097 reads processed maya sid all

# discarding tags that have more than 7 observations without reverse-complement
awk '!($3>7 && $4==0) && $2!="seq"' all.uniq >all.tab

wc -l all.tab  

# 147539 all.tab # before tufts custom code
# 177670 all.tab # after tufts custom code
# 230845 all.tab #maya
# 146282 all.tab #maya porites
# 77628 all.tab #maya sid
# 194224 all.tab maya sid combined
# 303541 all.tab maya sid all


# creating fasta file out of merged and filtered tags:
awk '{print ">"$1"\n"$2}' all.tab > all.fasta

# clustering reads into loci using cd-hit
# clustering allowing for up to 3 mismatches (-c 0.91); the most abundant sequence becomes reference
# -aL : alignment coverage for the longer sequence (setting at 1 means the alignment must cover 100% of the sequence)
# -aS : alignment coverage for the shorter sequence (setting to 1 means the alignment must cover 100% of the sequence
# -c : sequence identity threshold, calculated as # of identical bases in alignment / full length of the shorter sequence
# -g : the program will cluster it into the most similar cluster that meets the threshold, rather than the first
# -M -T : memory limit and # of threads to use (0 sets no limits)

module load blast
module load cdhit
#module load cdhit/4.6.8 - maya couldn't get this specific version, hoping it works ok

cd-hit-est -i all.fasta -o cdh_alltags.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 0 -T 0  

#------------
# making fake reference genome (of 30 chromosomes) out of cd-hit cluster representatives
# need bowtie2, samtools and picard_tools for indexing

module load bowtie2
module load samtools
module load picard
module load perl

/proj/kdcastil/users/mayapow/2bRAD/2bRAD_denovo/concatFasta.pl fasta=cdh_alltags.fas num=30

#hannah concatenating 99753 records into 30 pseudo-chromosomes 
#hannah 3326 records per chromosome

#concatenating 133064 records into 30 pseudo-chromosomes  maya sid all
#4436 records per chromosome

#concatenating 105988 records into 30 pseudo-chromosomes maya sid combo hannah
#3533 records per chromosome

# maya porites concatenating 84945 records into 30 pseudo-chromosomes 
# maya porites 2832 records per chromosome
# maya sid 1898 records per chromosome

#maya sid all concatenating 141349 records into 30 pseudo-chromosomes 
#4712 records per chromosome

# formatting fake genome

export GENOME_FASTA=cdh_alltags_cc.fasta

bowtie2-build $GENOME_FASTA $GENOME_FASTA

samtools faidx $GENOME_FASTA


#==============
# Mapping reads to denovo reference and formatting bam files 

# for denovo: map with bowtie2 with default parameters
# this step will take a long time
>maps
for F in `ls *.fastq`; do
echo "bowtie2 --no-unal -x $GENOME_FASTA -U $F -S $F.sam">>maps
done

bash maps

ls *.fastq.sam > sams
cat sams | wc -l  # number should match number of.fastq files
# 48 - matches maya
# 24 - matches maya porites & sid
#74 matches maya sid combined
#85 with hannah data
#130 all combined

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
cat sams | perl -pe 's/(\S+)\.sam/samtools view -bS $1\.sam >$1\.unsorted\.bam && samtools sort $1\.unsorted\.bam -o $1\.bam && samtools index $1\.bam/' >s2b

bash s2b

rm -f *unsorted*
ls *bam | wc -l  # should be the same number as number of .fastq files
# 48! 
#24
#85
#130

# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis. Archive them.
cd ..
mkdir bam_files_por
mkdir bam_files_sid
mkdir bam_files_combo
mkdir bam_files_all
#scp -r /proj/kdcastil/users/mayapow/2bRAD/denovo_nosym/*bam /proj/kdcastil/users/mayapow/2bRAD/bam_files
scp -r /proj/kdcastil/users/mayapow/2bRAD/denovo_por/*bam /proj/kdcastil/users/mayapow/2bRAD/bam_files_por
scp -r /proj/kdcastil/users/mayapow/2bRAD/denovo_sid/*bam /proj/kdcastil/users/mayapow/2bRAD/bam_files_sid
scp -r /proj/kdcastil/users/mayapow/2bRAD/denovo_sid_combo/*bam /proj/kdcastil/users/mayapow/2bRAD/bam_files_combo
scp -r /proj/kdcastil/users/mayapow/2bRAD/denovo_nosym/*bam /proj/kdcastil/users/mayapow/2bRAD/bam_files_all

##for moving bam files to local computer to use in r
tar -czvf bam_files_all.tar.gz /proj/kdcastil/users/mayapow/2bRAD/bam_files_all/*.bam
tar -xvf bam_files_all.tar.gz


#------------ quality assessment
###THIS DOES NOT WORK BC I DIDNT SUBMIT TO CLUSTER SO IT DOESNT MAKE maps.e from call above
#>alignmentRates
#>align
#for F in `ls *nosymbio.fastq`; do 
#echo "grep -E '^[ATGCN]+$' $F | wc -l | grep -f - maps.e* -A 4 | tail -1 | perl -pe 's/maps\.e\d+-|% overall alignment rate//g' | xargs echo $F.sam >> alignmentRates" >>align; 
#done
#bash align
# this is still not working for me...alignmentRates is just a list of the .sam file names after the job finishes running
# MAYA HAD THIS SAME ISSUE
#grep: maps.e*: No such file or directory

#FIXED!
#instead run this:
#$F.bam is the name your s2b step produced from $F.sam.
#bc we mapped with --no-unal, unaligned reads aren't in the BAM
#dividing mapped reads by the fastq read count therefore gives the same number as bowtie2's "overall alignment rate"
#-F 260 drops unmapped and secondary alignments.

>alignmentRates
for F in *.fastq; do
  total=$(( $(wc -l < $F) / 4 ))
  mapped=$(samtools view -c -F 260 $F.bam)
  awk -v f=$F.sam -v m=$mapped -v t=$total 'BEGIN{printf "%s\t%.2f\n", f, 100*m/t}' >> alignmentRates
done

cat alignmentRates

#por by itself
#18_S18_L001_R1_001_ACCA.nosymbio.fastq.sam      83.28
#18_S18_L001_R1_001_AGAC.nosymbio.fastq.sam      85.84
#18_S18_L001_R1_001_AGTG.nosymbio.fastq.sam      85.36
#18_S18_L001_R1_001_CATC.nosymbio.fastq.sam      83.59
#18_S18_L001_R1_001_CTAC.nosymbio.fastq.sam      86.40
#18_S18_L001_R1_001_GACT.nosymbio.fastq.sam      86.11
#18_S18_L001_R1_001_GCTT.nosymbio.fastq.sam      83.57
#18_S18_L001_R1_001_GTGA.nosymbio.fastq.sam      82.47
#18_S18_L001_R1_001_TCAC.nosymbio.fastq.sam      63.24
#18_S18_L001_R1_001_TCAG.nosymbio.fastq.sam      86.47
#18_S18_L001_R1_001_TGGT.nosymbio.fastq.sam      85.14
#18_S18_L001_R1_001_TGTC.nosymbio.fastq.sam      84.75
#19_S19_L001_R1_001_ACCA.nosymbio.fastq.sam      87.66
#19_S19_L001_R1_001_AGAC.nosymbio.fastq.sam      88.40
#19_S19_L001_R1_001_AGTG.nosymbio.fastq.sam      79.67
#19_S19_L001_R1_001_CATC.nosymbio.fastq.sam      77.07
#19_S19_L001_R1_001_CTAC.nosymbio.fastq.sam      87.98
#19_S19_L001_R1_001_GACT.nosymbio.fastq.sam      87.99
#19_S19_L001_R1_001_GCTT.nosymbio.fastq.sam      87.32
#19_S19_L001_R1_001_GTGA.nosymbio.fastq.sam      86.50
#19_S19_L001_R1_001_TCAC.nosymbio.fastq.sam      83.74
#19_S19_L001_R1_001_TCAG.nosymbio.fastq.sam      88.36
#19_S19_L001_R1_001_TGGT.nosymbio.fastq.sam      87.91
#19_S19_L001_R1_001_TGTC.nosymbio.fastq.sam      86.38

#sid by itself
#16_S16_L001_R1_001_ACCA.nosymbio.fastq.sam      77.14
#16_S16_L001_R1_001_AGAC.nosymbio.fastq.sam      76.53
#16_S16_L001_R1_001_AGTG.nosymbio.fastq.sam      77.45
#16_S16_L001_R1_001_CATC.nosymbio.fastq.sam      77.14
#16_S16_L001_R1_001_CTAC.nosymbio.fastq.sam      75.99
#16_S16_L001_R1_001_GACT.nosymbio.fastq.sam      75.59
#16_S16_L001_R1_001_GCTT.nosymbio.fastq.sam      75.02
#16_S16_L001_R1_001_GTGA.nosymbio.fastq.sam      76.71
#16_S16_L001_R1_001_TCAC.nosymbio.fastq.sam      75.38
#16_S16_L001_R1_001_TCAG.nosymbio.fastq.sam      76.81
#16_S16_L001_R1_001_TGGT.nosymbio.fastq.sam      76.55
#16_S16_L001_R1_001_TGTC.nosymbio.fastq.sam      75.79
#17_S17_L001_R1_001_ACCA.nosymbio.fastq.sam      69.64
#17_S17_L001_R1_001_AGAC.nosymbio.fastq.sam      69.67
#17_S17_L001_R1_001_AGTG.nosymbio.fastq.sam      71.07
#17_S17_L001_R1_001_CATC.nosymbio.fastq.sam      70.48
#17_S17_L001_R1_001_CTAC.nosymbio.fastq.sam      70.26
#17_S17_L001_R1_001_GACT.nosymbio.fastq.sam      69.61
#17_S17_L001_R1_001_GCTT.nosymbio.fastq.sam      69.44
#17_S17_L001_R1_001_GTGA.nosymbio.fastq.sam      68.99
#17_S17_L001_R1_001_TCAC.nosymbio.fastq.sam      68.29
#17_S17_L001_R1_001_TCAG.nosymbio.fastq.sam      70.38
#17_S17_L001_R1_001_TGGT.nosymbio.fastq.sam      66.80
#17_S17_L001_R1_001_TGTC.nosymbio.fastq.sam      69.80

#maya combo with hannah data
16_S16_L001_R1_001_ACCA.nosymbio.fastq.sam      85.24
16_S16_L001_R1_001_AGAC.nosymbio.fastq.sam      85.31
16_S16_L001_R1_001_AGTG.nosymbio.fastq.sam      86.07
16_S16_L001_R1_001_CATC.nosymbio.fastq.sam      85.32
16_S16_L001_R1_001_CTAC.nosymbio.fastq.sam      84.40
16_S16_L001_R1_001_GACT.nosymbio.fastq.sam      83.99
16_S16_L001_R1_001_GCTT.nosymbio.fastq.sam      83.50
16_S16_L001_R1_001_GTGA.nosymbio.fastq.sam      85.29
16_S16_L001_R1_001_TCAC.nosymbio.fastq.sam      84.00
16_S16_L001_R1_001_TCAG.nosymbio.fastq.sam      85.21
16_S16_L001_R1_001_TGGT.nosymbio.fastq.sam      84.78
16_S16_L001_R1_001_TGTC.nosymbio.fastq.sam      84.55
17_S17_L001_R1_001_ACCA.nosymbio.fastq.sam      80.79
17_S17_L001_R1_001_AGAC.nosymbio.fastq.sam      80.34
17_S17_L001_R1_001_AGTG.nosymbio.fastq.sam      82.08
17_S17_L001_R1_001_CATC.nosymbio.fastq.sam      81.30
17_S17_L001_R1_001_CTAC.nosymbio.fastq.sam      81.71
17_S17_L001_R1_001_GACT.nosymbio.fastq.sam      80.70
17_S17_L001_R1_001_GCTT.nosymbio.fastq.sam      80.39
17_S17_L001_R1_001_GTGA.nosymbio.fastq.sam      80.92
17_S17_L001_R1_001_TCAC.nosymbio.fastq.sam      80.39
17_S17_L001_R1_001_TCAG.nosymbio.fastq.sam      81.30
17_S17_L001_R1_001_TGGT.nosymbio.fastq.sam      81.96
17_S17_L001_R1_001_TGTC.nosymbio.fastq.sam      81.03
SRR24593955.fastq.sam   99.63
SRR24593956.fastq.sam   99.60
SRR24593957.fastq.sam   99.51
SRR24593958.fastq.sam   99.58
SRR24593959.fastq.sam   99.60
SRR24593960.fastq.sam   99.62
SRR24593961.fastq.sam   99.58
SRR24593962.fastq.sam   99.61
SRR24593963.fastq.sam   99.64
SRR24593964.fastq.sam   99.63
SRR24593965.fastq.sam   99.63
SRR24593966.fastq.sam   99.61
SRR24593967.fastq.sam   99.64
SRR24593968.fastq.sam   99.60
SRR24593969.fastq.sam   99.61
SRR24593970.fastq.sam   99.61
SRR24593971.fastq.sam   99.60
SRR24593972.fastq.sam   99.59
SRR24593973.fastq.sam   99.59
SRR24593974.fastq.sam   99.62
SRR24593975.fastq.sam   99.61
SRR24593976.fastq.sam   99.58
SRR24593977.fastq.sam   99.66
SRR24593978.fastq.sam   99.53
SRR24593979.fastq.sam   99.59
SRR24593980.fastq.sam   99.61
SRR24593981.fastq.sam   99.61
SRR24593982.fastq.sam   99.62
SRR24593983.fastq.sam   99.67
SRR24593984.fastq.sam   99.56
SRR24593985.fastq.sam   99.64
SRR24593986.fastq.sam   99.61
SRR24593987.fastq.sam   99.63
SRR24593988.fastq.sam   99.61
SRR24593989.fastq.sam   99.57
SRR24593990.fastq.sam   99.58
SRR24593991.fastq.sam   99.56
SRR24593992.fastq.sam   99.65
SRR24593993.fastq.sam   99.40
SRR24593994.fastq.sam   99.54
SRR24593995.fastq.sam   99.56
SRR24593996.fastq.sam   99.62
SRR24593997.fastq.sam   99.61
SRR24593998.fastq.sam   99.63
SRR24593999.fastq.sam   99.61
SRR24594000.fastq.sam   99.60
SRR24594001.fastq.sam   99.66
SRR24594002.fastq.sam   99.61
SRR24594003.fastq.sam   99.58
SRR24594004.fastq.sam   99.68
SRR24594005.fastq.sam   99.62
SRR24594006.fastq.sam   99.59
SRR24594007.fastq.sam   99.61
SRR24594008.fastq.sam   99.63
SRR24594009.fastq.sam   99.66
SRR24594010.fastq.sam   99.62
SRR24594011.fastq.sam   99.67
SRR24594012.fastq.sam   99.63
SRR24594013.fastq.sam   99.62
SRR24594014.fastq.sam   99.67
SRR24594015.fastq.sam   99.59

#------------------------------ANGSD (fuzzy genotyping)
# this first run-through of ANGSD is just to take a look at base qualities and coverage depth, 
# we will run angsd again with filters informed by this first run-through

pwd
#/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym
/proj/kdcastil/users/mayapow/2bRAD/denovo_por

# "FUZZY genotyping" with ANGSD - without calling actual genotypes but working with genotype likelihoods at each SNP. Optimal for low-coverage data (<10x).

# listing all bam filenames 
ls *bam >bams

#module load angsd
#used this link to install angsd using the unix instructions - put in my 2bRAD folder
#https://popgen.dk/angsd/index.php/Installation
#useable at: 
/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd


#don't have module on cluster so have to use from this folder instead
#make everything executable (navigate to angsd folder)
chmod +x *.cpp
chmod +x *.h
chmod +x *.c

#----------- assessing base qualities and coverage depth

# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping =< 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd : the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.

#FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 610 -minInd 30"
#FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 480 -minInd 24"
#FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 240 -minInd 12"
#FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 740 -minInd 37"
#FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 850 -minInd 42"
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1300 -minInd 65"
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

# in the following line, -r argument is one chromosome or contig to work with 
#(no need to do this for whole genome as long as the chosen chromosome or contig is long enough, ~1 Mb. 
#When mapping to a real genome, consider chr1:1-1000000 )
# (look up lengths of your contigs in the header of *.sam files)
# i'm not using here and just running on entire reference

#angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out dd
/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out dd


# summarizing results (using modified script by Matteo Fumagalli)
module load r
Rscript /proj/kdcastil/users/mayapow/2bRAD/2bRAD_denovo/plotQC.R prefix=dd

# proportion of sites covered at >5x:
cat quality.txt

#Porites - Maya - these look really good
#18_S18_L001_R1_001_TCAC.nosymbio.fastq.bam 0.0207417585087513
#18_S18_L001_R1_001_TGTC.nosymbio.fastq.bam 0.420398467954187
#18_S18_L001_R1_001_GACT.nosymbio.fastq.bam 0.540354097901664
#18_S18_L001_R1_001_CTAC.nosymbio.fastq.bam 0.576173239440074
#18_S18_L001_R1_001_CATC.nosymbio.fastq.bam 0.703392063040719
#18_S18_L001_R1_001_AGTG.nosymbio.fastq.bam 0.720991344627503
#19_S19_L001_R1_001_CATC.nosymbio.fastq.bam 0.727968780329547
#18_S18_L001_R1_001_ACCA.nosymbio.fastq.bam 0.741884452048491
#18_S18_L001_R1_001_GTGA.nosymbio.fastq.bam 0.744863973152743
#18_S18_L001_R1_001_TCAG.nosymbio.fastq.bam 0.777342943393391
#19_S19_L001_R1_001_GCTT.nosymbio.fastq.bam 0.793937507955331
#19_S19_L001_R1_001_GTGA.nosymbio.fastq.bam 0.81728609429724
#19_S19_L001_R1_001_CTAC.nosymbio.fastq.bam 0.83548420617785
#19_S19_L001_R1_001_TGTC.nosymbio.fastq.bam 0.838618107370362
#18_S18_L001_R1_001_GCTT.nosymbio.fastq.bam 0.84168577418216
#19_S19_L001_R1_001_TCAG.nosymbio.fastq.bam 0.850459715808889
#19_S19_L001_R1_001_GACT.nosymbio.fastq.bam 0.850944093135525
#18_S18_L001_R1_001_AGAC.nosymbio.fastq.bam 0.854640245439671
#19_S19_L001_R1_001_AGTG.nosymbio.fastq.bam 0.858678191576432
#19_S19_L001_R1_001_ACCA.nosymbio.fastq.bam 0.86182359832394
#19_S19_L001_R1_001_AGAC.nosymbio.fastq.bam 0.879159793711342
#19_S19_L001_R1_001_TCAC.nosymbio.fastq.bam 0.901034708808619
#18_S18_L001_R1_001_TGGT.nosymbio.fastq.bam 0.907167945640532
#19_S19_L001_R1_001_TGGT.nosymbio.fastq.bam 0.910579465908173

#Siderastrea - Maya - these look pretty bad for 1/2 of them (pretty much all bay I think :()
#16_S16_L001_R1_001_CTAC.nosymbio.fastq.bam      0.01632342
#16_S16_L001_R1_001_CATC.nosymbio.fastq.bam      0.01694571
#16_S16_L001_R1_001_ACCA.nosymbio.fastq.bam      0.01895616
#16_S16_L001_R1_001_GACT.nosymbio.fastq.bam      0.01981820
#16_S16_L001_R1_001_GCTT.nosymbio.fastq.bam      0.02079583
#16_S16_L001_R1_001_TGTC.nosymbio.fastq.bam      0.02204000
#16_S16_L001_R1_001_AGTG.nosymbio.fastq.bam      0.02266822
#16_S16_L001_R1_001_TGGT.nosymbio.fastq.bam      0.02290353
#16_S16_L001_R1_001_AGAC.nosymbio.fastq.bam      0.02380320
#16_S16_L001_R1_001_TCAG.nosymbio.fastq.bam      0.02586068
#16_S16_L001_R1_001_GTGA.nosymbio.fastq.bam      0.02671344
#16_S16_L001_R1_001_TCAC.nosymbio.fastq.bam      0.02717982
#17_S17_L001_R1_001_AGAC.nosymbio.fastq.bam      0.27013559
#17_S17_L001_R1_001_TCAG.nosymbio.fastq.bam      0.27071546
#17_S17_L001_R1_001_CATC.nosymbio.fastq.bam      0.30362765
#17_S17_L001_R1_001_CTAC.nosymbio.fastq.bam      0.37984905
#17_S17_L001_R1_001_AGTG.nosymbio.fastq.bam      0.43461092
#17_S17_L001_R1_001_TGGT.nosymbio.fastq.bam      0.46506344
#17_S17_L001_R1_001_TGTC.nosymbio.fastq.bam      0.47116390
#17_S17_L001_R1_001_ACCA.nosymbio.fastq.bam      0.53969602
#17_S17_L001_R1_001_GACT.nosymbio.fastq.bam      0.57822589
#17_S17_L001_R1_001_GCTT.nosymbio.fastq.bam      0.64692112
#17_S17_L001_R1_001_TCAC.nosymbio.fastq.bam      0.73164628
#17_S17_L001_R1_001_GTGA.nosymbio.fastq.bam      0.75404553

#Siderastrea with hannah's data
16_S16_L001_R1_001_CTAC.nosymbio.fastq.bam 0.0150236838240874
16_S16_L001_R1_001_CATC.nosymbio.fastq.bam 0.0173340097508504
16_S16_L001_R1_001_ACCA.nosymbio.fastq.bam 0.0183087550368356
16_S16_L001_R1_001_GACT.nosymbio.fastq.bam 0.0183541873267191
16_S16_L001_R1_001_GCTT.nosymbio.fastq.bam 0.0187098983112021
16_S16_L001_R1_001_AGTG.nosymbio.fastq.bam 0.0210057051043303
16_S16_L001_R1_001_TGTC.nosymbio.fastq.bam 0.0213303645170727
16_S16_L001_R1_001_TGGT.nosymbio.fastq.bam 0.0218489157047855
16_S16_L001_R1_001_AGAC.nosymbio.fastq.bam 0.0219466935331651
16_S16_L001_R1_001_TCAG.nosymbio.fastq.bam 0.0242828191237022
16_S16_L001_R1_001_GTGA.nosymbio.fastq.bam 0.0259501783645963
16_S16_L001_R1_001_TCAC.nosymbio.fastq.bam 0.027147798008334
SRR24593957.fastq.bam 0.0843693453115205
SRR24593984.fastq.bam 0.0991325532137045
SRR24593993.fastq.bam 0.110037682509703
SRR24593969.fastq.bam 0.137863259874324
SRR24593992.fastq.bam 0.139313805688615
SRR24593982.fastq.bam 0.147493576238866
SRR24594015.fastq.bam 0.149577889030772
SRR24594014.fastq.bam 0.152775234593817
SRR24593973.fastq.bam 0.161930337236636
SRR24593972.fastq.bam 0.162068460163254
SRR24594008.fastq.bam 0.167364332592173
SRR24593987.fastq.bam 0.170643678945448
SRR24593978.fastq.bam 0.172507783052496
SRR24594001.fastq.bam 0.173938758003079
SRR24593994.fastq.bam 0.176684285296348
SRR24593968.fastq.bam 0.17887832112506
SRR24593979.fastq.bam 0.182965958114069
SRR24594013.fastq.bam 0.191251175117913
SRR24593980.fastq.bam 0.192196465269456
SRR24594007.fastq.bam 0.196772612015003
SRR24593959.fastq.bam 0.206310328446489
SRR24594011.fastq.bam 0.209147026750225
SRR24593955.fastq.bam 0.214459438773864
SRR24593985.fastq.bam 0.220337780121009
SRR24593963.fastq.bam 0.220489965777582
SRR24593995.fastq.bam 0.22080466220609
SRR24593983.fastq.bam 0.246943526432452
SRR24593974.fastq.bam 0.25602775553034
SRR24593976.fastq.bam 0.261489441863257
SRR24593977.fastq.bam 0.262668698982588
SRR24594004.fastq.bam 0.266514204415172
17_S17_L001_R1_001_AGAC.nosymbio.fastq.bam 0.267354840348225
SRR24594006.fastq.bam 0.26877218787513
17_S17_L001_R1_001_TCAG.nosymbio.fastq.bam 0.269448083991644
17_S17_L001_R1_001_CATC.nosymbio.fastq.bam 0.300663972425571
SRR24593970.fastq.bam 0.305407878631017
SRR24593996.fastq.bam 0.33178687416317
SRR24594010.fastq.bam 0.332991575575585
SRR24593960.fastq.bam 0.334978723259608
SRR24593966.fastq.bam 0.340468389974794
SRR24594002.fastq.bam 0.346347925803765
SRR24594009.fastq.bam 0.357764806927858
SRR24593965.fastq.bam 0.358514741821104
SRR24593967.fastq.bam 0.36959351548043
17_S17_L001_R1_001_CTAC.nosymbio.fastq.bam 0.381420810033132
SRR24593964.fastq.bam 0.386542493511528
SRR24593989.fastq.bam 0.394014300673818
SRR24594005.fastq.bam 0.395131764741644
SRR24594003.fastq.bam 0.406212626740542
SRR24593997.fastq.bam 0.407518497478093
SRR24593971.fastq.bam 0.419269047503771
SRR24593998.fastq.bam 0.423988774993965
SRR24593975.fastq.bam 0.427718742640665
17_S17_L001_R1_001_AGTG.nosymbio.fastq.bam 0.440549378160845
17_S17_L001_R1_001_TGTC.nosymbio.fastq.bam 0.476907091187675
17_S17_L001_R1_001_TGGT.nosymbio.fastq.bam 0.48711504417937
SRR24593962.fastq.bam 0.494255516003677
SRR24593986.fastq.bam 0.5159321924913
SRR24593961.fastq.bam 0.523542280219748
17_S17_L001_R1_001_ACCA.nosymbio.fastq.bam 0.53672417421804
SRR24593988.fastq.bam 0.546340929217762
17_S17_L001_R1_001_GACT.nosymbio.fastq.bam 0.574703305804786
SRR24593991.fastq.bam 0.588200606988091
SRR24593958.fastq.bam 0.590318139808865
SRR24594000.fastq.bam 0.608125298870823
SRR24594012.fastq.bam 0.633036686626497
SRR24593956.fastq.bam 0.633086667467062
SRR24593999.fastq.bam 0.635590731021121
17_S17_L001_R1_001_GCTT.nosymbio.fastq.bam 0.643608378176063
17_S17_L001_R1_001_TCAC.nosymbio.fastq.bam 0.732923773604801
SRR24593990.fastq.bam 0.739010381638517
17_S17_L001_R1_001_GTGA.nosymbio.fastq.bam 0.756400047942418
SRR24593981.fastq.bam 0.838455322750344

# scp dd.pdf to laptop to look at distribution of base quality scores, 
#fraction of sites in each sample passing coverage thresholds, 
#and fraction of sites passing genotyping rates cutoffs. 
#Use these to guide choices of -minQ,  -minIndDepth and -minInd filters in subsequent ANGSD runs

# try again without the -baq, since this was a new filter for me - 
#JP Rippe used it but Misha doesn't have it in his run through
#mp got exact same results here whether you have -baq filter or not - no need!!

#----------- clones detection (**minInd ~80% of samples)

pwd
#/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym
/proj/kdcastil/users/mayapow/2bRAD/denovo_por

#FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 49 -snp_pval 1e-5 -minMaf 0.05"
#TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2"
#TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2 "
#SET minInd to 80% of samples!!
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 19 -snp_pval 1e-5 -minMaf 0.05"
TODO="-domajorminor 1 -domaf 1 -docounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -dogeno 8 -dobcf 1 -doPost 1 -doGlf 2"
#-doVcf is deprecated. Please use -doBcf 1 
# example : -dobcf 1 -doMajorMinor 1 -doPost 1 -gl 1 -domaf 1 -docounts 1 -dogeno 1 <- got this note while running, edited above accordingly
#[bcfoutput]     Please add the following parameters '-gl 1 -dopost 1 -domajorminor 1 -domaf 1 -dobcf 1 --ignore-RG 0 -dogeno 1 -docounts 1'

#for sid - bc such low quality, maya trying to lower minInd to see if technical replicates will show as clones
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 15 -snp_pval 1e-5 -minMaf 0.05"
TODO="-domajorminor 1 -domaf 1 -docounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -dogeno 8 -dobcf 1 -doPost 1 -doGlf 2"

#for sid combo
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 60 -snp_pval 1e-5 -minMaf 0.05"
TODO="-domajorminor 1 -domaf 1 -docounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -dogeno 8 -dobcf 1 -doPost 1 -doGlf 2"

#for sid all
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 95 -snp_pval 1e-5 -minMaf 0.05"
TODO="-domajorminor 1 -domaf 1 -docounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -dogeno 8 -dobcf 1 -doPost 1 -doGlf 2"


/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out myresult

NSITES=`zcat myresult.mafs.gz | wc -l`
echo $NSITES
# 11635 before tufts pipeline
# 14319 after tufts custom pipeline
#maya porites 16780
#maya sid 408 EEEEK with minInd = 19
#with minInd = 15 -> 6534 much better!!
#with minInd = 12 -> 25870 also good but want to stick closer to the 80% so using 15 since it still separates the clones
#919 with minInd = 59 - try lower maybe 50 sid combo
#11029 with minInd = 50 sid combo
#15332 with minInd = 60 sid hannah combo
#10584 with minInd = 95 sid all

# the -doVcf flag produces a vcf file at this stage.
# with this new version of angsd it doesn't work to create the vcf file - need to run:
bcftools view myresult.bcf -O z -o myresult.vcf.gz #this creates vcf file
bcftools view myresult.bcf -O v -o myresult.vcf #creates the uncompressed vcf file
# can also use: -O v if you want the uncompressed file to use 
#uncompressed is the file plink --vcf the ngsRelate and bayescan steps expect
tabix -p vcf myresult.vcf.gz #this creates the .tbi

# Add the -minIndDepth filter to see what happens, this filter was really important in Nicola's project but not included in Misha's pipeline. 
# Starting with -minIndDepth 2 and see what happens when we increase it
#  -minIndDepth is unknown will exit Maya - not using

# scp *.ibsMat and bams to laptop, use ibs.R to analyze

# On my local machine, the denovo without sym contamination files are located here:

pwd
/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/old_denovo_nosyms

# the denovo without sym contamination files plus the tufts custom pipeline files are here:

pwd
/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms


#------------------------------RE-RUN WITHOUT CLONES - Lineage assignments 

pwd
TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

# first with all duplicate preps removed, but with the two clones detected still included
# note the -minInd changes, to ~80% of our now 51 samples
# taking out all clones!

#files to remove that are clones or tech reps:
#chose tech reps with lower quality ofc
16_S16_L001_R1_001_TCAG.nosymbio.fastq.bam
16_S16_L001_R1_001_TGGT.nosymbio.fastq.bam
17_S17_L001_R1_001_GACT.nosymbio.fastq.bam
17_S17_L001_R1_001_TGTC.nosymbio.fastq.bam
SRR24593957.fastq.bam
SRR24593958.fastq.bam
SRR24593960.fastq.bam
SRR24593961.fastq.bam
SRR24593962.fastq.bam
SRR24593963.fastq.bam
SRR24593964.fastq.bam
SRR24593966.fastq.bam
SRR24593975.fastq.bam
SRR24593976.fastq.bam
SRR24593979.fastq.bam

#should be 70 left (-15)
#put into text file : remove.txt and uploaded to cluster
#had issues with leftover characters in remove.txt - edited accordingly with this:
sed -i 's/\r$//; s/[[:space:]]*$//' remove.txt


grep -v -x -F -f remove.txt bams > bams_noclones

wc -l bams_noclones

#now with clones removed - only 113 samples left so 

#SET minInd to 80% of samples!!
#por
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 18 -snp_pval 1e-5 -minMaf 0.05"
TODO="-domajorminor 1 -domaf 1 -docounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -dogeno 8 -dobcf 1 -doPost 1 -doGlf 2"

#sid
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 80 -snp_pval 1e-5 -minMaf 0.05"
TODO="-domajorminor 1 -domaf 1 -docounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -dogeno 8 -dobcf 1 -doPost 1 -doGlf 2"

/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd -b bams_noclones -GL 1 $FILTERS $TODO -P 1 -out myresult.noclones

NSITES=`zcat myresult.noclones.mafs.gz | wc -l`
echo $NSITES
# 7176
# 14878 MP hannah combine
# 14716 MP hannah combine noclones
# 12947 MP sid all w/min ind 80 for 113 samples
# 15646 MP por noclones

# the -doVcf flag produces a vcf file at this stage.
# with this new version of angsd it doesn't work to create the vcf file - need to run:
bcftools view myresult.noclones.bcf -O z -o myresult.noclones.vcf.gz #this creates vcf file
bcftools view myresult.noclones.bcf -O v -o myresult.noclones.vcf #creates the uncompressed vcf file
# can also use: -O v if you want the uncompressed file to use 
#uncompressed is the file plink --vcf the ngsRelate and bayescan steps expect
tabix -p vcf myresult.noclones.vcf.gz #this creates the .tbi

#download myresult.noclones.ibsMat

#------------------------------NgsAdmix

# first, install:
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/bin

wget popgen.dk/software/download/NGSadmix/ngsadmix32.cpp 

# compile
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix

# all set, downloaded in bin.
#Useable for maya at:
/proj/kdcastil/users/mayapow/2bRAD/NGSadmix

# NgsAdmix for K from 2 to 5 : do not run if the dataset contains clones or genotyping replicates!

#for K in `seq 2 5` ;
#do
#../../bin/NGSadmix -likes myresult2.noclone.beagle.gz -K $K -P 10 -o mydata.noclone_k${K};
#done

for K in `seq 1 5` ;
do
/proj/kdcastil/users/mayapow/2bRAD/NGSadmix -likes myresult.noclones.beagle.gz -K $K -P 10 -o mydata_noclones_k${K};
done

#scp *qopt files to local machine, use admixturePlotting_v5.R to plot


##also use real ADMIXTURE on called SNPs (requires plink and ADMIXTURE):

pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

module load admixture/1.3.0
#module load plink/1.90b6.4 - couldnt get this version
module load plink
# unzip vcf file:
#gunzip myresult.vcf.gz
#gunzip myresult.vcf.gz
module list
#check versions using above

#already unzipped myresult.vcf - use here

# remove the "(angsd version)" from the head of the file, this is what was causing errors before

#plink --vcf myresult2.noclone.vcf --make-bed --allow-extra-chr 0 --out myresult2.noclone

#plink was for old model now plink2
#plink2 --vcf myresult.vcf --make-bed --chr-set 95 no-xy no-mt --allow-extra-chr 0 --double-id --out myresult
plink2 --vcf myresult.noclones.vcf --make-bed --chr-set 95 no-xy no-mt --allow-extra-chr 0 --double-id --out myresult.noclones
#issue with extra chromosome still

#check that files were made
#ls myresult.bed myresult.bim myresult.fam

#and that bim file looks as expected
#head myresult.bim

#for K in `seq 1 5`; \
#do admixture --cv myresult2.noclone.bed $K | tee myresult2.noclone_${K}.out; done


#plink --vcf myresult2.noclone.vcf --make-bed --allow-extra-chr 0 --out myresult2.noclone

#with clones
#for K in `seq 1 5`; \
#do admixture --cv myresult.bed $K | tee myresult_admix${K}.out; done

#and no clones
for K in `seq 1 5`; \
do admixture --cv myresult.noclones.bed $K | tee myresult.noclones_admix${K}.out; done

#Use this to check K of least CV error:
#grep -h "CV error" myresult_admix*.out
grep -h "CV error" myresult.noclones_admix*.out

#maya por
#CV error (K=1): 0.60727
#CV error (K=2): 0.43157
#CV error (K=3): 0.34114
#CV error (K=4): 0.36052
#CV error (K=5): 0.46157

#maya sid hannah
#CV error (K=1): 0.50310
#CV error (K=2): 0.44073
#CV error (K=3): 0.45133
#CV error (K=4): 0.47057
#CV error (K=5): 0.50333

#maya sid hannah no clones
#CV error (K=1): 0.50445
#CV error (K=2): 0.45494
#CV error (K=3): 0.47940
#CV error (K=4): 0.52524
#CV error (K=5): 0.60331

#seems like either way k = 2 or k = 3 is best

#maya sid all
#CV error (K=1): 0.51004
#CV error (K=2): 0.43361
#CV error (K=3): 0.43744
#CV error (K=4): 0.45047
#CV error (K=5): 0.47139

#------------------------------Analysis of Genetic Divergence

pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

# with all clones removed, we have 50 samples

#module load angsd/0.923
/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd

# filtering sites to work on - use only filters that do not distort allele frequency
# set minInd to 75-90% of the total number fo individuals in the project 
#maya set to 56 for 70 sample dataset - at 80%
# if you are doing any other RAD than 2bRAD or GBS, remove '-sb_pval 1e-5' from FILTERS
# -maxHetFreq 0.5 is the lumped paralog filter - important for doing SFS analyses
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 56"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doBcf 1"


#angsd -b bams_noclones -GL 1 $FILTERS $TODO -out sfilt
/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd -b bams_noclones -GL 1 $FILTERS $TODO -out sfilt

NSITES=`zcat sfilt.mafs.gz | wc -l`
echo $NSITES
#1967730 mp por
#393583 sid hannah no clones

# -maxHetFreq 0.5 is the lumped paralog filter - important for doing SFS analyses. seeing if it makes a difference for Fst values
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 56 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doBcf 1"


/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd -b bams_noclones -GL 1 $FILTERS $TODO -out sfilt_maxhet

NSITES=`zcat sfilt_maxhet.mafs.gz | wc -l`
echo $NSITES
# 770398
# 1964781 mp por
# 392994 sid hannah no clones

# make vcf only for running bayescan:
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 56 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doBcf 1"

/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd -b bams_noclones -GL 1 $FILTERS $TODO -out sfilt_forbayescan

NSITES=`zcat sfilt_forbayescan.mafs.gz | wc -l`
echo $NSITES
#15384 mp por
#4992 mp sid hannah no clones

# collecting and indexing filter-passing sites
zcat sfilt.mafs.gz | cut -f 1,2 | tail -n +2 >allSites
/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd sites index allSites

# max het filter
zcat sfilt_maxhet.mafs.gz | cut -f 1,2 | tail -n +2 >allSites_maxhet
/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd sites index allSites_maxhet

export GENOME_REF=/proj/kdcastil/users/mayapow/2bRAD/denovo_sid_combo/cdh_alltags_cc.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"

#need to do this by lineage now once I talk with sarah to discuss this!!
#MP TO DO:
#- go through pipeline and run with all belize samples and Hannah samples together
/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd -sites allSites -b bams_noclones -GL 1 -P 1 $TODO -out mapp
angsd -sites allSites -b bams_L2 -GL 1 -P 1 $TODO -out L2
angsd -sites allSites -b bams_L3 -GL 1 -P 1 $TODO -out L3

# max het filter
/proj/kdcastil/users/mayapow/2bRAD/angsd/angsd -sites allSites_maxhet -b bams -GL 1 -P 1 $TODO -out mapp_maxhet

angsd -sites allSites_maxhet -b bams_L1 -GL 1 -P 1 $TODO -out L1_maxhet
angsd -sites allSites_maxhet -b bams_L2 -GL 1 -P 1 $TODO -out L2_maxhet
angsd -sites allSites_maxhet -b bams_L3 -GL 1 -P 1 $TODO -out L3_maxhet

# generating per-population SFS
realSFS L1.saf.idx >L1.sfs
realSFS L2.saf.idx >L2.sfs
realSFS L3.saf.idx >L3.sfs

/proj/kdcastil/users/mayapow/2bRAD/angsd/realSFS mapp.saf.idx >mapp.sfs

# generating per-population SFS - max het filter
realSFS L1_maxhet.saf.idx >L1_maxhet.sfs
realSFS L2_maxhet.saf.idx >L2_maxhet.sfs
realSFS L3_maxhet.saf.idx >L3_maxhet.sfs

# writing down 2d-SFS priors - L1 vs L2
realSFS L1_maxhet.saf.idx L2_maxhet.saf.idx -P 24 > p12_maxhet.sfs ; realSFS fst index L1_maxhet.saf.idx L2_maxhet.saf.idx -sfs p12_maxhet.sfs -fstout p12_maxhet

# global Fst between populations
realSFS fst stats p12_maxhet.fst.idx

#output:
	-> Assuming idxname:p12.fst.idx
	-> Assuming .fst.gz file: p12.fst.gz
	-> FST.Unweight[nObs:771675]:0.023862 Fst.Weight:0.159775
0.023862 0.159775

	-> Assuming idxname:p12_maxhet.fst.idx
	-> Assuming .fst.gz file: p12_maxhet.fst.gz
	-> FST.Unweight[nObs:770220]:0.023784 Fst.Weight:0.172208
0.023784	0.172208

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].


# writing down 2d-SFS priors - L1 vs L3
realSFS L1_maxhet.saf.idx L3_maxhet.saf.idx -P 24 > p13_maxhet.sfs ; realSFS fst index L1_maxhet.saf.idx L3_maxhet.saf.idx -sfs p13_maxhet.sfs -fstout p13_maxhet

# global Fst between populations
realSFS fst stats p13_maxhet.fst.idx

#output:
	-> Assuming idxname:p13.fst.idx
	-> Assuming .fst.gz file: p13.fst.gz
	-> FST.Unweight[nObs:761904]:0.326915 Fst.Weight:0.162287
0.326915 0.162287

	-> Assuming idxname:p13_maxhet.fst.idx
	-> Assuming .fst.gz file: p13_maxhet.fst.gz
	-> FST.Unweight[nObs:760462]:0.327435 Fst.Weight:0.178868
0.327435	0.178868

# writing down 2d-SFS priors - L2 vs L3
realSFS L2_maxhet.saf.idx L3_maxhet.saf.idx -P 24 > p23_maxhet.sfs ; realSFS fst index L2_maxhet.saf.idx L3_maxhet.saf.idx -sfs p23_maxhet.sfs -fstout p23_maxhet

# global Fst between populations
realSFS fst stats p23_maxhet.fst.idx

#output:
	-> Assuming idxname:p23.fst.idx
	-> Assuming .fst.gz file: p23.fst.gz
	-> FST.Unweight[nObs:761878]:0.178944 Fst.Weight:0.108837
0.178944 0.108837

	-> Assuming idxname:p23_maxhet.fst.idx
	-> Assuming .fst.gz file: p23_maxhet.fst.gz
	-> FST.Unweight[nObs:760437]:0.178861 Fst.Weight:0.117511
0.178861	0.117511

# using the max het filter here (removing lumped paralogs) since it is most appropriate for SFS analyses

#------------------------------Identifying Outlier Loci
# To do this, will need to download PGDspider and Bayescan

#----- PGDspider :

cd ~/bin
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/bin

wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.0.7.1.zip
unzip PGDSpider_2.0.7.1.zip

#----- Bayescan :

cd ~/bin
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/bin

wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
unzip BayeScan2.1.zip
cp BayeScan2.1/binaries/BayeScan2.1_linux64bits bayescan
chmod +x bayescan
rm -r -f BayeScan*


#============= Bayescan: looking for Fst outliers
## AGAIN ISSUE HERE IS VCF FILE

# Converting vcf (using PGDspider) to Bayescan format: 
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

# make tab-delimited file called bspops LISTING assignments of individuals (as they are named in the vcf file) to populations, for example:
ind1	pop0
ind2	pop0
ind3	pop1
ind4	pop1

# create a file called vcf2bayescan.spid containing this text:
nano vcf2bayescan.spid

echo "############
# VCF Parser questions
PARSER_FORMAT=VCF
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=true
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=true
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=./bspops
# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=false
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# GESTE / BayeScan Writer questions
WRITER_FORMAT=GESTE_BAYE_SCAN
# Specify which data type should be included in the GESTE / BayeScan file  (GESTE / BayeScan can only analyze one data type per file):
GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP
############" >vcf2bayescan.spid

# launching bayescan (Misha says this might take 12-24 hours so submitting as a job)
head bayescan

#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N bayescan # job name, anything you want
#$ -l h_rt=48:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

java -Xmx49152m -Xms512m -jar /projectnb/davies-hb/hannah/TVE_2bRAD/bin/PGDSpider_2.0.7.1/PGDSpider2-cli.jar -inputfile sfilt_forbayescan.vcf -outputfile Best.bayescan -spid vcf2bayescan.spid 

/projectnb/davies-hb/hannah/TVE_2bRAD/bin/bayescan Best.bayescan -threads=20

# then submit job:
qsub -pe omp 28 bayescan


# use bayescan_plots.R to examine results

# removing outliers from VCF file in order to re-calculate Fst without outliers (neutral sites only)
# still have to use sfilt here!

#Best.baye_fst.txt is your bayescan output
cut -d" " -f2- Best.baye_fst.txt | tail -n +2> aaa
#find line number in vcf file that contains header
grep -n "#CHROM" sfilt_forbayescan.vcf
#13, add 1
tail --lines=+14 sfilt_forbayescan.vcf | cut -f 1,2 | paste --delimiters "\t" - aaa > baye_fst_pos.txt

awk '$5<0.5 {print $1"\t"$2}' ./baye_fst_pos.txt > bayesOut4FST

# bayeOuts4FST is the list of outlier loci
grep -Fwvf bayesOut4FST sfilt.vcf > NeutralLoci.vcf
grep -Fwf bayesOut4FST sfilt.vcf > OutlierLoci.vcf


zgrep "^[^#]" NeutralLoci.vcf | awk '{print $1,$2}' >NeutralLoci.sites
zgrep "^[^#]" OutlierLoci.vcf | awk '{print $1,$2}' >OutlierLoci.sites


#------------------------------Analysis of Genetic Divergence Again - Neutral Loci Only
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

# with all clones removed, we have 50 samples

module load angsd/0.923

# collecting and indexing filter-passing sites but without outliers
# cut -f 1,2 sfilt.nooutlier.vcf > Sites_noOutliers
# remove the header lines from the Sites_noOutliers file before indexing
angsd sites index NeutralLoci.sites

export GENOME_REF=/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym/cdh_alltags_cc.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"


angsd -sites NeutralLoci.sites -b bams_L1 -GL 1 -P 1 $TODO -out L1_neutral
angsd -sites NeutralLoci.sites -b bams_L2 -GL 1 -P 1 $TODO -out L2_neutral
angsd -sites NeutralLoci.sites -b bams_L3 -GL 1 -P 1 $TODO -out L3_neutral

# generating per-population SFS
realSFS L1_neutral.saf.idx >L1_neutral.sfs
realSFS L2_neutral.saf.idx >L2_neutral.sfs
realSFS L3_neutral.saf.idx >L3_neutral.sfs

# writing down 2d-SFS priors - L1 vs L2
realSFS L1_neutral.saf.idx L2_neutral.saf.idx -P 1 > p12_neutral.sfs ; realSFS fst index L1_neutral.saf.idx L2_neutral.saf.idx -sfs p12_neutral.sfs -fstout p12_neutral

# global Fst between populations
realSFS fst stats p12_neutral.fst.idx

#output:
	-> Assuming idxname:p12_neutral.fst.idx
	-> Assuming .fst.gz file: p12_neutral.fst.gz
	-> FST.Unweight[nObs:771324]:0.023573 Fst.Weight:0.136383
0.023573 0.136383

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].


# writing down 2d-SFS priors - L1 vs L3
realSFS L1_neutral.saf.idx L3_neutral.saf.idx -P 1 > p13_neutral.sfs ; realSFS fst index L1_neutral.saf.idx L3_neutral.saf.idx -sfs p13_neutral.sfs -fstout p13_neutral

# global Fst between populations
realSFS fst stats p13_neutral.fst.idx

#output:
	-> Assuming idxname:p13_neutral.fst.idx
	-> Assuming .fst.gz file: p13_neutral.fst.gz
	-> FST.Unweight[nObs:761562]:0.325370 Fst.Weight:0.143440
0.325370 0.143440

# writing down 2d-SFS priors - L2 vs L3
realSFS L2_neutral.saf.idx L3_neutral.saf.idx -P 1 > p23_neutral.sfs ; realSFS fst index L2_neutral.saf.idx L3_neutral.saf.idx -sfs p23_neutral.sfs -fstout p23_neutral

# global Fst between populations
realSFS fst stats p23_neutral.fst.idx

#output:
	-> Assuming idxname:p23_neutral.fst.idx
	-> Assuming .fst.gz file: p23_neutral.fst.gz
	-> FST.Unweight[nObs:761536]:0.178447 Fst.Weight:0.100000
0.178447 0.100000


#------------------------------Analysis of Genetic Divergence Again - Outlier Loci Only
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

# with all clones removed, we have 50 samples

module load angsd/0.923

# remove the header lines from the Sites_noOutliers file before indexing
angsd sites index OutlierLoci.sites

export GENOME_REF=/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym/cdh_alltags_cc.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"


angsd -sites OutlierLoci.sites -b bams_L1 -GL 1 -P 1 $TODO -out L1_outliers
angsd -sites OutlierLoci.sites -b bams_L2 -GL 1 -P 1 $TODO -out L2_outliers
angsd -sites OutlierLoci.sites -b bams_L3 -GL 1 -P 1 $TODO -out L3_outliers

# generating per-population SFS
realSFS L1_outliers.saf.idx >L1_outliers.sfs
realSFS L2_outliers.saf.idx >L2_outliers.sfs
realSFS L3_outliers.saf.idx >L3_outliers.sfs

# writing down 2d-SFS priors - L1 vs L2
realSFS L1_outliers.saf.idx L2_outliers.saf.idx -P 1 > p12_outliers.sfs ; realSFS fst index L1_outliers.saf.idx L2_outliers.saf.idx -sfs p12_outliers.sfs -fstout p12_outliers

# global Fst between populations
realSFS fst stats p12_outliers.fst.idx

#output:
	-> Assuming idxname:p12_outliers.fst.idx
	-> Assuming .fst.gz file: p12_outliers.fst.gz
	-> FST.Unweight[nObs:351]:0.645093 Fst.Weight:0.732279
0.645093 0.732279

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].


# writing down 2d-SFS priors - L1 vs L3
realSFS L1_outliers.saf.idx L3_outliers.saf.idx -P 1 > p13_outliers.sfs ; realSFS fst index L1_outliers.saf.idx L3_outliers.saf.idx -sfs p13_outliers.sfs -fstout p13_outliers

# global Fst between populations
realSFS fst stats p13_outliers.fst.idx

#output:
	-> Assuming idxname:p13_outliers.fst.idx
	-> Assuming .fst.gz file: p13_outliers.fst.gz
	-> FST.Unweight[nObs:342]:0.659837 Fst.Weight:0.713172
0.659837 0.713172

# writing down 2d-SFS priors - L2 vs L3
realSFS L2_outliers.saf.idx L3_outliers.saf.idx -P 1 > p23_outliers.sfs ; realSFS fst index L2_outliers.saf.idx L3_outliers.saf.idx -sfs p23_outliers.sfs -fstout p23_outliers

# global Fst between populations
realSFS fst stats p23_outliers.fst.idx

#output:
	-> Assuming idxname:p23_outliers.fst.idx
	-> Assuming .fst.gz file: p23_outliers.fst.gz
	-> FST.Unweight[nObs:342]:0.466759 Fst.Weight:0.603811
0.466759 0.603811


#------------------------------DOWNSTREAM ANALYSES

# Used angsd_ibs_pca.R to plot hierarchical clustering of samples based on identity by state matrix (file = myresult2.noclone.ibsMat)
# in addition to principal coordinate analysis (PCoA) also using the same ibsMat.
# Used admixturePlotting_v5.R to plot results of NgsAdmix (file = mydata.noclone_k3.qopt)


#------------------------------NGSRelate
# Running NGSRelate during revision to check on relatedness of L3 individuals:
# Additional workflow details here: https://github.com/ANGSD/NgsRelate
# Download and Install

pwd
#/projectnb/davies-hb/hannah/TVE_Panama/TVE_2bRAD/bin
/proj/kdcastil/users/mayapow/2bRAD

git clone --recursive https://github.com/SAMtools/htslib
git clone https://github.com/ANGSD/ngsRelate
cd htslib/;make -j2;cd ../ngsRelate;make HTSSRC=../htslib/
#make everything executable (navigate to ngsRelate folder)
#chmod +x *.cpp
#chmod +x *.h

#useable at:
/proj/kdcastil/users/mayapow/2bRAD/ngsRelate/ngsRelate

pwd
#/projectnb/davies-hb/hannah/TVE_Panama/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym
/proj/kdcastil/users/mayapow/2bRAD/denovo_nosym

#gzipped vcf file with samples excluding technical replicates and clones:
myresult.noclone.vcf.gz
myresult.vcf.gz

#/projectnb/davies-hb/hannah/TVE_Panama/TVE_2bRAD/bin/ngsRelate/ngsRelate -h /projectnb/davies-hb/hannah/TVE_Panama/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym/myresult2.noclone.vcf.gz.tbi -O ngsrelate_vcf.res

/proj/kdcastil/users/mayapow/2bRAD/ngsRelate/ngsRelate -h myresult.vcf.gz -O ngsrelate_vcf.res

# doing this on two datasets - one with technical replicates and clones removed (ngsrelate.noclone) and one with technical replicates removed but one clone pair retained (ngsrelate.noclone.allsamps)
# following the pipeline on the ngsrelate website: https://github.com/ANGSD/NgsRelate

#MAYA NOT DOING ANY OF THIS BELOW BC GOT ABOVE TO WORK

### First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).
angsd -b bams_noclones -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -out ngsrelate.noclone
NSITES=`zcat ngsrelate.noclone.mafs.gz | wc -l`
echo $NSITES
# 104277

angsd -b bams_noclones_allsamps -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -out ngsrelate.noclone.allsamps
NSITES=`zcat ngsrelate.noclone.allsamps.mafs.gz | wc -l`
echo $NSITES
# 104205

### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat ngsrelate.noclone.mafs.gz | cut -f5 |sed 1d >freq.noclone

zcat ngsrelate.noclone.allsamps.mafs.gz | cut -f5 |sed 1d >freq.noclone.allsamps

### run NgsRelate
/projectnb/davies-hb/hannah/TVE_Panama/TVE_2bRAD/bin/ngsRelate/ngsRelate -g ngsrelate.noclone.glf.gz -n 50 -f freq.noclone -O ngsrelate.noclone.res

/projectnb/davies-hb/hannah/TVE_Panama/TVE_2bRAD/bin/ngsRelate/ngsRelate -g ngsrelate.noclone.allsamps.glf.gz -n 51 -f freq.noclone.allsamps -O ngsrelate.noclone.allsamps.res

# copy *.res output to local computer. output of ngsrelate is analyzed in the angsd_ibs_pca.R script
