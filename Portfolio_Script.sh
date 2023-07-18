#! /bin/bash

chmod a+x Portfolio_Script.sh

#######################################################
########### PREAMBLE FOR ALL SHELL SCRIPTS ############
#######################################################
### Author: T. Williams								###
### Bioinformatics Pipeline							###
### Last Updated: July 18, 2023						###
#######################################################
#######################################################
### This file contains the tools and commands used in the entire variant calling pipline taking a raw .fastq file from sequencing and taking it all the way through picard 
### This file contains the steps starting with preparing the directory for processing and ends with Picard.
### This script is the first script that should be run. 
### Upon completion of this script, execute XXXX.sh for phase two of Variant analysis.
#######################################################
### Steps Overview in this pipeline
### Step 1: Define path variables for tool packages
### Step 2A: Define current working directory and create path variables for input and output directories
### Step 2B: Create an array for use in loop actions that creates your input and output directories
### Step 2C: Define path variable shortcuts for input and output directories
### Step 2D: Move raw .fastq files from main directory into starting directory
### Step 3A: use fastqc to assess quality on individual raw .fastq files
### Step 3B: use multiqc to assess quality on all raw .fastq files

### NEED TO DO THE FOLLOWING THINGS:
### RUN ALL FILES ALL THE WAY THROUGH TO BAM USING ZIPPED FILES. NO NEED TO UNZIP ANYTHING.

### PIPE THE OUTPUT OF BWA INTO SAMTOOLS SORT AND GET A SORTED BAM OUTPUT TO SAVE TIME


#######################################################
################# Begin File Comments #################
#######################################################
### Hello, welcome to the USDA Wooden Breast Variant Calling pipeline, USDA WB Quality Assessment Preprocessing Version 1.0. script. 
### 
### This is the first script that should be run in this GATK variant analysis pipeline suggested by the Broad Institute as this completes the entire data preprocessing pipeline.
### 
### The starting file is a raw .fastq files from sequencing and produces the necessary QC documents and statistics for whole genome sequenicng."
### 
### GATK workflow following completed using this script 
### 
### Set up preamble and transfer 28 practice .fastq files. 
### 
#######################################################
################## End File Comments ##################
#######################################################

#******# Step 1: Define path variables #******#
#******# Establish your current working directory #******#
cwdir=/Travis
echo ""
echo "Current working directory has been established as: $cwdir"
echo ""

#******# Define path variables for server wide tool packages #******#
#pkgs=/home/texasam/Bioinformatics/pkgs
picard=$cwdir/pkgs/picard.jar
#gatk=$pkgs/GATK405/gatk
#snpeff=$pkgs/snpEff/snpEff.jar
#qualimap=$pkgs/qualimap/qualimap
#bedtools=$pkgs/bedtools2

#******# Specify path variable for reference assembly #******#
gallus7=$cwdir/RefGenomes/Galgal7.fa

echo ""
echo "the reference genome is set to:" $gallus7
echo ""


#******# Create an array for use in loop actions that creates your input and output directories #******#
declare -a arr=(
	"00_SnL" 
	"01_RawFastq" 
	"02_Rawqc" 
	"03_Merged"
	"04_Trimmed" 
	"05_Posttrim_qc"
	"06_Sam" 
	"07_Bam" 
	"08_Alignment_Stats" 
	"09_Realigned" 
	"10_Rawvcf" 
	"11_Stats"
	"12_Snps" 
	"13_SnpEff"
	"14_Delly"
	"15_Mappability"
	"16_Additional3"
	"17_Additional4"
	"18_Additional5"
	"19_Additional6"
	"20_DataDownloads"
	)


#******# Define path variable shortcuts for input and output directories #******#
SnL=$cwdir/00_SnL && echo "Path variable set for 00_SnL:" $SnL 
RawFastq=$cwdir/01_RawFastq && echo "Path variable set for 01_RawFastq:" $RawFastq 
Rawqc=$cwdir/02_Rawqc && echo "Path variable set for 02_Rawqc:" $Rawqc 
Merged=$cwdir/03_Merged && echo "Path variable set for 03_Merged:" $Merged
Trimmed=$cwdir/04_Trimmed && echo "Path variable set for 04_Trimmed:" $Trimmed 
Posttrim_qc=$cwdir/05_Posttrim_qc && echo "Path variable set for 05_Posttrim_qc:" $Posttrim_qc 
Sam=$cwdir/06_Sam && echo "Path variable set for 06_Sam:" $Sam 
Bam=$cwdir/07_Bam && echo "Path variable set for 07_Bam:" $Sam 
Alignment_Stats=$cwdir/08_Alignment_Stats && echo "Path variable set for 08_Alignment_Stats:" $Alignment_Stats 
Realigned=$cwdir/09_Realigned && echo "Path variable set for 09_Realigned:" $Realigned 
Rawvcf=$cwdir/10_Rawvcf && echo "Path variable set for 10_Rawvcf:" $Rawvcf 
VCFStats=$cwdir/11_Stats && echo "Path variable set for 11_Stats:" $VCFStats
Snps=$cwdir/12_Snps && echo "Path variable set for 12_Snps:" $Snps 
SnpEff=$cwdir/13_SnpEff && echo "Path variable set for 13_SnpEff:" $SnpEff
Delly=$cwdir/14_Delly && echo "Path variable set for 14_Delly:" $Delly
Map=$cwdir/15_Mappability && echo "Path variable set for 15_Mappability:" $Map
Additional3=$cwdir/16_Additional3 && echo "Path variable set for 16_Additional3:" $Additional3
Additional4=$cwdir/17_Additional4 && echo "Path variable set for 17_Additional4:" $Additional4
Additional5=$cwdir/18_Additional5 && echo "Path variable set for 18_Additional5:" $Additional5
Additional6=$cwdir/19_Additional6 && echo "Path variable set for 19_Additional6:" $Additional6
DataDownloads=$cwdir/20_DataDownloads && echo "Path variable set for 20_DataDownloads:" $DataDownloads


#######################################################
#################### END PREAMBLE #####################
#######################################################


#######################################################
############### BEGIN SCRIPT COMMANDS #################
#######################################################

quit ( ) { # quit function that quits the script
	exit
}
## Update
echo ""
echo "Hello, welcome to the Portfolio Script bioinformatics pipeline."
echo "FYI, the sample file name and starting step needs to be passed in this format: bash Portfolio_Script.sh file.txt answer"
echo "This script can perform these steps by using them in the answer section of running this script:"
echo "DirSetup = Set up the directory structure to perform this analysis. This only needs to be run once at the beginning of the analysis."
echo "DirReset = Reset the directory structure. This is if you ever need to start over."
echo "Rename = Rename files to standardize nomenclature for latter use."
echo "RawQC = Performs QC on raw, untrimmed fastq files."
echo "FCLanes = Generates a list of flowcell lanes to verify sample consolidating in the next step. "
echo "CombineRawFQ = Combines rawfasq files for further use. "
echo "TrimReport = generate a report of trimming action. "
echo "RefIndex = Index the reference genome. "
echo "RapidAlign = rapidly aligns trimmed files to the reference genome using bwa-mem2. "
echo "QMReportConsolidation = generates a single file report of important data that is mined from qualimap reports. "
echo "PicardMD = Marks duplicates using the picard program. "
echo "DellyPrep = step to confirm that all necessary files needed to run delly have been generated. "
echo "DellyCall = calls structural variants using the delly program. "
echo "vcftrim = trims vcf files and cleans up data. "
echo "VariantCalling = calls snp variants"
echo "SnpFiltering = Filters snps called. "
echo "Exit = Exit this script."
echo "Type <begin text> quit <end text> at anytime to immediately quits this script."

samplelist=$1
answer=$2

echo 
echo "Your file containing sample names [blank if not given]:" $samplelist 
echo ""
echo "Your option to pass into this script along with your sample file [blank if not given]:" $answer
echo ""

#******# Begin script steps to run #******#

if [[ $answer == DirSetup ]]  
	then
		echo "you chose $answer, Creating Directory Structure"
			#******# Run once when setting up the master directory #******#
			for i in "${arr[@]}"
				do 
					mkdir -p $cwdir/$i
			done

			echo ""
			echo "Input and output directories have been created."
			echo ""

		echo ""	
		echo "Step $answer completed. Now Exiting..."
		echo ""

elif [[ $answer == DirReset ]] 
	then
		echo "you chose $answer, Resetting Directory Structure"		
			#******# Only run this step if you need to remove the directories you just created #******# 
			for i in "${arr[@]}"
				do 
					rm -rf $i
			done

			echo ""
			echo "Input and output directories have been removed."
			echo ""

		echo "Step $answer completed. Exiting..."
			
elif [[ $answer == Rename ]] 
	then
		echo ""
		echo "you chose $answer, beginning the combination step of raw fq.gz files on a per sample, per read, per lane basis..."
		echo ""
			#******# Comments:
			#******# Original file names are in the format of FlowcellID_Breed-Gen-Band-PM_SampleNumber_LaneNumber_ReadNumber_001.fq.gz 
			#******# Regardless of which run number it is. which is a problem as there is 2 flowcells, Flowcell1 and Flowcell2
			#******#
			#******# New File names need to be in Breed_Gen_Band_PM_FlowcellID_LaneNumber_ReadNumber.fq.gz
			#******#
			#******# run this script in a directory that has your fastq files for one flow cell and a "mapping"
			#******# file named "Flowcell1SampleNames.txt". "samplenames.txt" needs to be a tab delimited text file with
			#******# no headers. Importantly, the old file names need to be in column 1, and the new file names
			#******# need to be in column 2

			# Breed3ck if a filename for the mapping file is provided
			if [ $# -ne 1 ]; then
    			echo "Usage: $0 FlowcellSampleNames.txt"
    			exit 1
			fi

			# Read the mapping file line by line
			while IFS=$'\t' read -r oldname newname
			do
			    # Breed3ck if the file exists and is a regular file
			    if [ -f "$oldname" ]; then
			        # Rename the file
			        mv "$oldname" "$newname"
			        echo "Renamed $oldname to $newname"
			    else
			        echo "File $oldname not found!"
			    fi
			done < "$1"

		echo ""
		echo "Step $answer completed on all input files. Now Exiting..."
		echo ""			

elif [[ $answer == RawQC ]] 
	then
		echo ""
		echo "you chose $answer, beginning QC analysis on raw .fq files direct from the sequencing facility.."
		echo ""

		#******# use fastq on gzipped .fq files direct from sequencing facility. #******#
		while read sample;
			do 
			echo ""
			echo "Beginning FastQC analysis on file $sample"
			echo ""
			fastqc -t 6 $RawFastq/Flowcell1/$sample --outdir $Rawqc/Flowcell1/
			echo ""
			echo "Done performing FastQC on file $sample. Report has been written to $Rawqc/Flowcell1/"
			echo ""
		done < $samplelist

		for file in $RawFastq/Flowcell2/*.gz
			do 
			echo ""
			echo "Beginning FastQC analysis on $file"
			echo ""
			fastqc -t 6 $file --outdir $Rawqc/Flowcell2/
			echo ""
			echo "Done performing FastQC on $file. Report has been written to $Rawqc/Flowcell2/"
			echo ""
		done

		echo ""
		echo "FastQC Reports for all raw .fq files have been written to:" $Rawqc
		echo ""

		#******# Use multiqc to assess quality on raw .fastq files on a per sample basis #******# 
		while read sample; do
			echo ""
			echo "Beginning MultiQC analysis on all raw .fastq files for $sample"
			echo ""
			multiqc $Rawqc/${sample%}*.fq.gz --filename ${sample%}_multiqc_report --outdir $Rawqc
			echo ""
			echo "Done with MultiQC analysis on all raw .fastq files for $sample"
			echo ""
		done < $samplelist
		echo ""
		echo "Done performing MultiQC analysis on all raw fastq files on a per sample basis."
		echo ""

elif [[ $answer == FCLanes ]]
	then
		echo ""
		echo "you chose $answer, beginning..."
		echo ""
		#******# This lists all of the R1 variables, separates out all of the different lanes, sorts them; 
		#******# then prints out the unique lists into the files R1_FC1_Lanes, R2_FC1_Lanes, R1_FC2_Lanes , R2_FC2_Lanes;
		#******# This information is used in the CombineRawFQ step. 
		ls -1 $RawFastq/Flowcell1/*.fastq.gz | awk -F '_' 'OFS="_" {print $7, $8, $9}' | sort | uniq | awk 'NR <=8'> $SnL/Flowcell1_Lanes.txt
		ls -1 $RawFastq/Flowcell2/*.fastq.gz | awk -F '_' 'OFS="_" {print $7, $8, $9}' | sort | uniq | awk 'NR <=8'> $SnL/Flowcell2_Lanes.txt
		echo ""
		echo "Step $answer completed. Now Exiting..."
		echo ""

elif [[ $answer == CombineRawFQ ]] 
	then
		echo ""
		echo "you chose $answer, beginning the combination step of raw fq.gz files on a per sample, per read, per lane basis..."
		echo ""

			while read sample; do
				#******# Combine all R1 lanes for flowcells Flowcell1 and Flowcell2 #******#
				echo ""
				echo "Combine all R1 lanes for $sample from all flowcells"
				echo ""
				cat \
				$RawFastq/Flowcell1/${sample%}_PM_N23059_L001_R1.fastq.gz \
				$RawFastq/Flowcell1/${sample%}_PM_N23059_L002_R1.fastq.gz \
				$RawFastq/Flowcell1/${sample%}_PM_N23059_L003_R1.fastq.gz\
				$RawFastq/Flowcell1/${sample%}_PM_N23059_L004_R1.fastq.gz \
				$RawFastq/Flowcell2/${sample%}_PM_N23060_L001_R1.fastq.gz \
				$RawFastq/Flowcell2/${sample%}_PM_N23060_L002_R1.fastq.gz \
				$RawFastq/Flowcell2/${sample%}_PM_N23060_L003_R1.fastq.gz\
				$RawFastq/Flowcell2/${sample%}_PM_N23060_L004_R1.fastq.gz > $Merged/${sample%}_R1.fastq.gz
				echo ""
				echo "Beginning FastQC analysis on $Merged/${sample%}_R1.fastq.gz"
				echo ""
				fastqc -t 10 $Merged/${sample%}_R1.fastq.gz --outdir $Merged/QC/
				echo ""
				echo "Done running FastQC on ${sample%}_R1.fastq.gz"
				echo ""
				#******# Combine all R2 lanes for flowcells Flowcell1 and Flowcell2 #******#
				echo ""
				echo "Combine all R2 lanes for $sample from all flowcells"
				echo ""
				cat \
				$RawFastq/Flowcell1/${sample%}_PM_N23059_L001_R2.fastq.gz \
				$RawFastq/Flowcell1/${sample%}_PM_N23059_L002_R2.fastq.gz \
				$RawFastq/Flowcell1/${sample%}_PM_N23059_L003_R2.fastq.gz \
				$RawFastq/Flowcell1/${sample%}_PM_N23059_L004_R2.fastq.gz \
				$RawFastq/Flowcell2/${sample%}_PM_N23060_L001_R2.fastq.gz \
				$RawFastq/Flowcell2/${sample%}_PM_N23060_L002_R2.fastq.gz \
				$RawFastq/Flowcell2/${sample%}_PM_N23060_L003_R2.fastq.gz \
				$RawFastq/Flowcell2/${sample%}_PM_N23060_L004_R2.fastq.gz > $Merged/${sample%}_R2.fastq.gz
				echo ""
				echo "Beginning FastQC analysis on $Merged/${sample%}_R2.fastq.gz"
				echo ""
				fastqc -t 10 $Merged/${sample%}_R2.fastq.gz --outdir $Merged/QC/
				echo ""
				echo "Done running FastQC on ${sample%}_R2.fastq.gz"
				echo ""
				echo ""
				echo "Done combining all fastq.gz files for $sample and running fastqc"
				echo ""
			done < $samplelist

		echo ""
		echo "Step $answer completed on all input files. Now Exiting..."
		echo ""

#******# Trimming step removed following sample consolidation to preserve IP #******#

elif [[ $answer == TrimReport ]] 
	then
		echo ""
		echo "you chose $answer, consolidating information from all of the trimming reports"
		echo ""
			for trimreport in $Trimmed/*trimming_report.txt; do
				awk 'NR==4 { fielda = $3 }
			     NR==26 { fieldb =$4 }
			     NR==27 { fieldc = $4; fieldd = $5 }
			     NR==28 { fielde = $5; fieldf = $6 }
			     NR==30 { fieldg = $4 }
			     NR==31 { fieldh = $2; fieldi = $4 }
			     NR==32 { fieldj = $4; fieldk = $6 }
			     END { print fielda "\t" fieldb "\t" fieldc "\t" fieldd "\t" fielde "\t" fieldf "\t" fieldg "\t" fieldh "\t" fieldi "\t" fieldj "\t" fieldk }' $trimreport
			done

		echo ""
		echo "Step $answer completed. Now Exiting..."
		echo ""	

elif [[ $answer == RefIndex ]]
	then
		echo ""
		echo "you chose $answer, now starting the indexing of the $gallus7 genome"
		echo ""
		bwa-mem2 index $gallus7
		echo ""
		echo "Done with indexing the gallus7 reference genome using bwa-mem2. Now exiting..."
		echo ""	

elif [[ $answer == RapidAlign ]] 
	then
		echo ""
		echo "you chose the step $answer, now starting alignment of trimmed.fq.gz files to the $gallus6 genome"
		echo ""
		while read sample; do 
			#******#  Begin Comments: #******# 
			#******# Perform alignments with bwa on filtered data    
			#******# The alignments take the format of:
			#******# bwa mem RefGenome SampleID_ReadF1.fq.gz SampleID_ReadF2.fq.gz > SampleID.sam
			#******# in this case, we are piping straight into a sorted bam file
			#******#  End Comments: #******# 
			#******#  Align the file #******# 
			echo ""
			echo "Beginning alignment of $sample at $(date +%x_%r)"
			bwa-mem2 mem -t 10 $gallus7 $Trimmed/${sample%}_R1_trimmed.fq.gz $Trimmed/${sample%}_R2_trimmed.fq.gz | samtools sort -@10 -T $Bam/${sample%}_tempfile -O bam -o $Bam/${sample%}_sorted.bam
			echo ""
			echo "Alignment, conversion, and sorting of ${sample%}_sorted.bam was completed at $(date +%x_%r)."
			echo ""
		done < $samplelist

		while read sample; do 
			#******# Index the .bam file #******#
			echo "Beginning Indexing on ${sample%}_sorted.bam at $(date +%x_%r)"
			samtools index $Bam/${sample%}_sorted.bam
			echo "${sample%}_sorted.bam has been sorted and indexed at $(date +%x_%r)"
			#******# Quality analysis of newly generated file #******#
			echo "Beginning samtools flagstat on ${sample%}_sorted.bam at $(date +%x_%r)."
			samtools flagstat $Bam/${sample%}_sorted.bam -O tsv
			echo "Beginning Qualimap quality analysis of ${sample%}_sorted.bam at $(date +%x_%r)."
			mkdir $Alignment_Stats/${sample%}_QM_Report/
			#******# qualimap bamqc is a very ram intensive process. Always Breed3ck ram allotment using the free command #******#
			qualimap bamqc -bam $Bam/${sample%}_sorted.bam -nt 10 --java-mem-size=32G -outdir $Alignment_Stats/${sample%}_QM_Report/
			echo "Quality has been assessed via flagstat and qualimap. Done with $sample at $(date +%x_%r). Moving on to the next one..."
		done < $samplelist

		echo ""
		echo "Step $answer completed on all input files. Now Exiting..."
		echo ""

elif [[ $answer == QMReportConsolidation ]]
	then
		echo ""
		echo "you chose $answer, beginning..."
		echo ""

		# Define the path to the output file
		output_file=$SnL/QMReportSummary.txt
	
		# Extract the lines and save them to the output file as columns
		echo "Bam File,Number of Reads,Mapped Reads,Percent Mapped Reads,Mean Mapping Quality,Error Rate,Number of MismatBreed3s,Number of Insertions,Insertion Percentage,Number of Deletions,Deletion Percentage,Homopolymer Indels,Mean Coverage,STD Coverage" > "$output_file"	

		while IFS= read -r sample; do 
			# Iterate over the 10 subdirectories
  			for ReportDir in $Alignment_Stats/${sample%}_QM_Report; do
			    # Set the path to the genome_results.txt file for the current sample and subdirectory
			    genome_file="$ReportDir/genome_results.txt"

			    # Get the various elements from different lines using awk and put the output into a .txt file that can be converted to CSV for summary analysis
			    awk 'NR == 6 { Bamfile = $4 }
					 NR == 20 { NoReads = $5 }
				     NR == 21 { MappedReads = $6; PCTMapped = $7 }
				     NR == 46 { MeanMapQual = $5 }
				     NR == 62 { ErrorRate = $5 }
				     NR == 63 { MismatBreed3s = $5 }
				     NR == 64 { Insertions = $5 }
				     NR == 65 { PCTInsertion = $7 }
				     NR == 66 { Deletions = $5 }
				     NR == 67 { PCTDeletion = $7 }
				     NR == 68 { HomopolymerIndels = $4 }
				     NR == 73 { Coverage = $4 }
				     NR == 74 { STDCoverage = $4 }
				     END { print Bamfile "\t" NoReads "\t" MappedReads "\t" PCTMapped "\t" MeanMapQual "\t" ErrorRate "\t" MismatBreed3s "\t" Insertions "\t" PCTInsertion "\t" Deletions "\t" PCTDeletion "\t" HomopolymerIndels "\t" Coverage "\t" STDCoverage}' "$genome_file" >> "$output_file"
			done
	  	done < $samplelist

		echo ""
		echo "Step $answer completed on all input files in $samplelist. Now Exiting..."
		echo ""

elif [[ $answer == PicardMD ]]
	then
		echo ""
		echo "you chose step $answer, beginning Picard Mark Duplicates..."
		echo ""
		while read sample; do
			#******# Step 11A: Mark duplicates in BAM files using picard # make sure your .bam file is indexed.
			echo "Sample: $Bam/${sample%}_sorted.bam found at $(date +%x_%r). Marking duplicates now..."
			java -jar $picard MarkDuplicates \
				--INPUT $Bam/${sample%}_sorted.bam \
				--OUTPUT $Bam/${sample%}_md.bam \
				--METRICS_FILE $Bam/${sample%}_metrix.txt \
				--USE_JDK_INFLATER true \
				--USE_JDK_DEFLATER true
			echo "Mark duplicates on ${sample%}_sorted.bam completed at $(date +%x_%r). Moving onto the next sample..."
			done < $samplelist
		echo ""
		echo "Picard Mark Duplicates completed. All samples processed in the file $samplelist are ready for MDindex step."
		echo ""
		echo ""
		echo "Step $answer completed on all input files in $samplelist. Now Exiting..."
		echo ""

#******# Three additional Picard steps removed following PicardMD to preserve IP #******#

elif [[ $answer == DellyPrep ]]  
	then
		echo ""
		echo "you chose step $answer. To run the DELLY program, it requires a sorted, indexed and duplicate marked"
		echo "bam file for every input sample."
		echo "Additionally, an indexed reference genome is required to identify split-reads."
		echo "Beginning DELLY Prep analysis now..."
		echo ""
		# Define the output file for missing files
		Missing_Delly_Files="$SnL/Missing_Delly_Files.txt"

		# Define the output file for ready files
		Ready_Delly_Files="$SnL/Ready_Delly_Files.txt"

		# Breed3ck if each file is present for each sample in the list
		while IFS= read -r sample; do
		    sorted_bamfile="_sorted.bam"
		    indexed_sorted_bamfile="_sorted.bam.bai"
		    md_bamfile="_md.bam"
		    indexed_md_bamfile="_md.bam.bai"
		    md_srg_bamfile="_md_SRG.bam"
		    indexed_md_srg_bamfile="_md_SRG.bai"

		    missing_files=()

		    if [[ ! -f "$Bam/bamfiles/${sample%}$sorted_bamfile" ]]; then
		        missing_files+=("${sample%}$sorted_bamfile")
		    fi

		    if [[ ! -f "$Bam/bamfiles/${sample%}$indexed_sorted_bamfile" ]]; then
		        missing_files+=("${sample%}$indexed_sorted_bamfile")
		    fi

		    if [[ ! -f "$Bam/mdbamfiles/${sample%}$md_bamfile" ]]; then
		        missing_files+=("${sample%}$md_bamfile")
		    fi

		    if [[ ! -f "$Bam/mdbamfiles/${sample%}$indexed_md_bamfile" ]]; then
		        missing_files+=("${sample%}$indexed_md_bamfile")
		    fi

		    if [[ ! -f "$Bam/SRGbamfiles/${sample%}$md_srg_bamfile" ]]; then
		        missing_files+=("${sample%}$md_srg_bamfile")
		    fi

		    if [[ ! -f "$Bam/SRGbamfiles/${sample%}$indexed_md_srg_bamfile" ]]; then
		        missing_files+=("${sample%}$indexed_md_srg_bamfile")
		   	fi

			if [[ ${#missing_files[@]} -eq 0 ]]; then
		        echo "All files are present for sample: $sample"
		        echo "Sample $sample is ready for DELLY Analysis." >> "$Ready_Delly_Files"
		    else
		        echo "Missing files for sample: $sample"
		        echo "Missing files: ${missing_files[*]}"
		        printf "%s\n" "${missing_files[@]}" >> "$Missing_Delly_Files"
		    fi

		done < $samplelist
		echo ""
		echo "DELLY preparation analysis completed on all samples processed in the file $samplelist."
		echo "Samples ready for DELLY analysis can be found in $Ready_Delly_Files."
		echo "Samples NOT READY for DELLY analysis can be found in $Missing_Delly_Files."
		echo ""
		echo "Step $answer completed. Now Exiting..."
		echo ""	

elif [[ $answer == DellyCall ]] # Run per breed, uncomment based on which breed youre running. 
	then
		echo ""
		echo "you chose step $answer. To run the DELLY program, it requires a sorted, indexed and duplicate marked"
		echo "bam file for every input sample. An indexed reference genome is required to identify split-reads."
		echo "Delly primarily parallelizes on the sample level."
		echo "Hence, OMP_NUM_THREADS should be always smaller or equal to the number of input samples."
		echo "Beginning DELLY analysis now..."
		echo ""
		# Delly supports parallel computing using the OpenMP API (www.openmp.org).
		#make PARALLEL=1 -B delly # This code does not apply to usage in a conda environment.
		#******# Just use the export OMP_NUM_THREADS code for each of the breeds.

		# You can set the number of threads using the environment variable OMP_NUM_THREADS.
		#******# Need to set this for each run.
		#******# for Breed1 use below.
		#export OMP_NUM_THREADS=35
		#******# For Breed3 use below.
		export OMP_NUM_THREADS=21
		#******# For Breed2 use below.
		#export OMP_NUM_THREADS=34

		while read sample; do
			#******# Delly call
			#******# -g indexed reference genome
			#******# input.bam is sample_md_SRG.bam file
			#******# What is -x?
			#echo "Beginning Delly analysis on ${sample%}_md_SRG.bam at $(date +%x_%r)."
			#******# for Breed1 use below.
			#delly call -g $gallus7 -o $Delly/Breed1/${sample%}.bcf $Bam/SRGbamfiles/${sample%}_md_SRG.bam &
			#******# Running 15:56 on 6.26.23
			#******# For Breed3 use below.
			delly call -g $gallus7 -o $Delly/Breed3/${sample%}.bcf $Bam/SRGbamfiles/${sample%}_md_SRG.bam &
			#******# Running  on 
			#******# For Breed2 use below.
			#delly call -g $gallus7 -o $Delly/Breed2/${sample%}.bcf $Bam/SRGbamfiles/${sample%}_md_SRG.bam &
			#******# Running 19:20 on 6.26.23
			#echo "Done with sample ${sample%}_md_SRG.bam. Moving onto the next sample."
		done < $samplelist
		#echo ""
		#echo "DELLY operations completed on all Samples. Next Steps are to merge, re-genotype, and filter the SV's..."
		#echo ""
		#echo ""
		#echo "Step $answer completed. Now Exiting..."
		#echo ""

#******# Additional DELLY steps removed to preserve IP #******#

elif [[ $answer == vcftrim ]] 
	then
		echo ""
		echo "you chose $answer, beginning the removal of unnecessary chromosomes at $(date +%x_%r)...."
		echo ""
		vcftools --vcf AS_germline.vcf --recode --recode-INFO-all --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr 23 --chr 24 --chr 25 --chr 26 --chr 27 --chr 28 --chr 29 --chr 30 --chr 31 --chr 32 --chr 33 --chr 34 --chr 35 --chr 36 --chr 37 --chr 38 --chr 39 --chr W --chr Z --chr MT  --remove-filtered-all --out AS_germline_trimmed
		
		#******# Generate Subsets of AS_germline_trimmed by breed
		vcftools --vcf AS_germline_trimmed.vcf --out Breed1_tgsv --recode --recode-INFO-all --indv Breed1_F1_001 --indv Breed1_F1_002 --indv Breed1_F1_004 --indv Breed1_F1_005 --indv Breed1_F1_020 --indv Breed1_F1_023 --indv Breed1_F1_025 --indv Breed1_F1_026 --indv Breed1_F1_028 --indv Breed1_F1_029 --indv Breed1_F1_031 --indv Breed1_F1_035 --indv Breed1_F1_038 --indv Breed1_F1_039 --indv Breed1_F1_040 --indv Breed1_F1_041 --indv Breed1_F1_050 --indv Breed1_F1_054 --indv Breed1_F1_060 --indv Breed1_F1_062 --indv Breed1_PS_H002 --indv Breed1_PS_H005 --indv Breed1_PS_H007 --indv Breed1_PS_H011 --indv Breed1_PS_H015 --indv Breed1_PS_H016 --indv Breed1_PS_H017 --indv Breed1_PS_H023 --indv Breed1_PS_H026 --indv Breed1_PS_H028 --indv Breed1_PS_R006 --indv Breed1_PS_R015 --indv Breed1_PS_R016 --indv Breed1_PS_R023 --indv Breed1_PS_RSTR_H012 
		vcftools --vcf AS_germline_trimmed.vcf --out Breed3_tgsv --recode --recode-INFO-all --indv Breed3_F1_001 --indv Breed3_F1_002 --indv Breed3_F1_003 --indv Breed3_F1_007 --indv Breed3_F1_008 --indv Breed3_F1_011 --indv Breed3_F1_012 --indv Breed3_F1_013 --indv Breed3_F1_015 --indv Breed3_F1_019 --indv Breed3_F1_022 --indv Breed3_F1_023 --indv Breed3_F1_024 --indv Breed3_F1_025 --indv Breed3_PS_H166 --indv Breed3_PS_H170 --indv Breed3_PS_H189 --indv Breed3_PS_H194 --indv Breed3_PS_R004 --indv Breed3_PS_R005 --indv Breed3_PS_R009 
		vcftools --vcf AS_germline_trimmed.vcf --out Breed2_tgsv --recode --recode-INFO-all --indv Breed2_F1_001 --indv Breed2_F1_002 --indv Breed2_F1_008 --indv Breed2_F1_012 --indv Breed2_F1_017 --indv Breed2_F1_018 --indv Breed2_F1_028 --indv Breed2_F1_032 --indv Breed2_F1_035 --indv Breed2_F1_063 --indv Breed2_F1_120 --indv Breed2_F1_128 --indv Breed2_F1_164 --indv Breed2_F1_166 --indv Breed2_F1_179 --indv Breed2_F1_183 --indv Breed2_F1_184 --indv Breed2_F1_185 --indv Breed2_F1_187 --indv Breed2_F1_195 --indv Breed2_PS_H084 --indv Breed2_PS_H088 --indv Breed2_PS_H098 --indv Breed2_PS_H343 --indv Breed2_PS_H345 --indv Breed2_PS_H355 --indv Breed2_PS_H357 --indv Breed2_PS_H368 --indv Breed2_PS_H373 --indv Breed2_PS_H383 --indv Breed2_PS_R001 --indv Breed2_PS_R002 --indv Breed2_PS_R003 --indv Breed2_PS_R007 
		
		#******# Following the delly germline filtering pipeline, remove all scaffolds and create a new trimmed germline vcf file by breed.
		vcftools --vcf Breed1_germline.vcf --recode --recode-INFO-all --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr 23 --chr 24 --chr 25 --chr 26 --chr 27 --chr 28 --chr 29 --chr 30 --chr 31 --chr 32 --chr 33 --chr 34 --chr 35 --chr 36 --chr 37 --chr 38 --chr 39 --chr W --chr Z --chr MT  --remove-filtered-all --out Breed1_germline_tsv
		vcftools --vcf Breed3_germline.vcf --recode --recode-INFO-all --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr 23 --chr 24 --chr 25 --chr 26 --chr 27 --chr 28 --chr 29 --chr 30 --chr 31 --chr 32 --chr 33 --chr 34 --chr 35 --chr 36 --chr 37 --chr 38 --chr 39 --chr W --chr Z --chr MT  --remove-filtered-all --out Breed3_germline_tsv
		vcftools --vcf Breed2_germline.vcf --recode --recode-INFO-all --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr 23 --chr 24 --chr 25 --chr 26 --chr 27 --chr 28 --chr 29 --chr 30 --chr 31 --chr 32 --chr 33 --chr 34 --chr 35 --chr 36 --chr 37 --chr 38 --chr 39 --chr W --chr Z --chr MT  --remove-filtered-all --out Breed2_germline_tsv

		
		#******# Subset by Generation
		vcftools --gzvcf Breed1_germline_tsv.recode.vcf.gz --out Breed1_germline_F1_tsv --recode --recode-INFO-all --indv Breed1_F1_001 --indv Breed1_F1_002 --indv Breed1_F1_004 --indv Breed1_F1_005 --indv Breed1_F1_020 --indv Breed1_F1_023 --indv Breed1_F1_025 --indv Breed1_F1_026 --indv Breed1_F1_028 --indv Breed1_F1_029 --indv Breed1_F1_031 --indv Breed1_F1_035 --indv Breed1_F1_038 --indv Breed1_F1_039 --indv Breed1_F1_040 --indv Breed1_F1_041 --indv Breed1_F1_050 --indv Breed1_F1_054 --indv Breed1_F1_060 --indv Breed1_F1_062 
		vcftools --gzvcf Breed1_germline_tsv.recode.vcf.gz --out Breed1_germline_PS_tsv --recode --recode-INFO-all --indv Breed1_PS_H002 --indv Breed1_PS_H005 --indv Breed1_PS_H007 --indv Breed1_PS_H011 --indv Breed1_PS_H015 --indv Breed1_PS_H016 --indv Breed1_PS_H017 --indv Breed1_PS_H023 --indv Breed1_PS_H026 --indv Breed1_PS_H028 --indv Breed1_PS_R006 --indv Breed1_PS_R015 --indv Breed1_PS_R016 --indv Breed1_PS_R023 --indv Breed1_PS_RSTR_H012 
		vcftools --gzvcf Breed3_germline_tsv.recode.vcf.gz --out Breed3_germline_F1_tsv --recode --recode-INFO-all --indv Breed3_F1_001 --indv Breed3_F1_002 --indv Breed3_F1_003 --indv Breed3_F1_007 --indv Breed3_F1_008 --indv Breed3_F1_011 --indv Breed3_F1_012 --indv Breed3_F1_013 --indv Breed3_F1_015 --indv Breed3_F1_019 --indv Breed3_F1_022 --indv Breed3_F1_023 --indv Breed3_F1_024 --indv Breed3_F1_025 
		vcftools --gzvcf Breed3_germline_tsv.recode.vcf.gz --out Breed3_germline_PS_tsv --recode --recode-INFO-all --indv Breed3_PS_H166 --indv Breed3_PS_H170 --indv Breed3_PS_H189 --indv Breed3_PS_H194 --indv Breed3_PS_R004 --indv Breed3_PS_R005 --indv Breed3_PS_R009 
		vcftools --gzvcf Breed2_germline_tsv.recode.vcf.gz --out Breed2_germline_F1_tsv --recode --recode-INFO-all --indv Breed2_F1_001 --indv Breed2_F1_002 --indv Breed2_F1_008 --indv Breed2_F1_012 --indv Breed2_F1_017 --indv Breed2_F1_018 --indv Breed2_F1_028 --indv Breed2_F1_032 --indv Breed2_F1_035 --indv Breed2_F1_063 --indv Breed2_F1_120 --indv Breed2_F1_128 --indv Breed2_F1_164 --indv Breed2_F1_166 --indv Breed2_F1_179 --indv Breed2_F1_183 --indv Breed2_F1_184 --indv Breed2_F1_185 --indv Breed2_F1_187 --indv Breed2_F1_195 
		vcftools --gzvcf Breed2_germline_tsv.recode.vcf.gz --out Breed2_germline_PS_tsv --recode --recode-INFO-all --indv Breed2_PS_H084 --indv Breed2_PS_H088 --indv Breed2_PS_H098 --indv Breed2_PS_H343 --indv Breed2_PS_H345 --indv Breed2_PS_H355 --indv Breed2_PS_H357 --indv Breed2_PS_H368 --indv Breed2_PS_H373 --indv Breed2_PS_H383 --indv Breed2_PS_R001 --indv Breed2_PS_R002 --indv Breed2_PS_R003 --indv Breed2_PS_R007 
		
		#******# Compress the VCF files for use with bcftools isec
		bgzip -@ 10 -k Breed2_germline_F1_WBD0_tsv.recode.vcf
		bgzip -@ 10 -k Breed2_germline_F1_WBD1_tsv.recode.vcf
		bgzip -@ 10 -k Breed2_germline_F1_WBD2_tsv.recode.vcf
		bgzip -@ 10 -k Breed2_germline_F1_WBD3_tsv.recode.vcf
		bgzip -@ 10 -k Breed3_germline_F1_WBD0_tsv.recode.vcf
		bgzip -@ 10 -k Breed3_germline_F1_WBD1_tsv.recode.vcf
		bgzip -@ 10 -k Breed3_germline_F1_WBD2_tsv.recode.vcf

		#******# Index the compressed vcf files for use with bcftools isec
		tabix Breed2_germline_F1_WBD0_tsv.recode.vcf.gz
		tabix Breed2_germline_F1_WBD1_tsv.recode.vcf.gz
		tabix Breed2_germline_F1_WBD2_tsv.recode.vcf.gz
		tabix Breed2_germline_F1_WBD3_tsv.recode.vcf.gz
		tabix Breed3_germline_F1_WBD0_tsv.recode.vcf.gz
		tabix Breed3_germline_F1_WBD1_tsv.recode.vcf.gz
		tabix Breed3_germline_F1_WBD2_tsv.recode.vcf.gz

		echo ""
		echo "Step $answer completed at $(date +%x_%r).. Now Exiting..."
		echo ""

elif [[ $answer == VariantCalling ]] 
	then
		echo ""
		echo "you chose $answer, beginning..."
		echo ""
		#******# Performing GATK Operations
		echo "Starting GATK Variant Calling steps"
		while read sample; do 
			echo "Sample" $sample "found"
			echo "Sample" $sample "InDel Realignment Target Creation"
			# d. Generate IndelRealigner Targets on Reference genome using GATK. NB: if you are working with a homogeneous sample set, only one interval file may suffice.  
			java -jar $gatk -T RealignerTargetCreator -R $gallus -I $bam/$sample"_SRG.bam" -o $bam/$sample".forIndelRealigner.intervals"
			echo "Sample" $sample "InDel Target Created"
			# e. Run IndelRealigner to realign original BAM files over INDEL sites using GATK
			java -jar $gatk -T IndelRealigner -R $gallus  -I $bam/$sample"_SRG.bam" --targetIntervals $bam/$sample".forIndelRealigner.intervals" -o $bam/$sample"_realigned.bam" --filter_bases_not_stored
			echo "Sample" $sample "realigned"
			# f. Run HaplotypeCaller on realigned BAM using GATK (variant discovery)
			java -jar $gatk -T HaplotypeCaller -R $gallus -I $bam/$sample"_realigned.bam" -o $rawvcf/$sample".raw.vcf" --min_base_quality_score 28 --maxNumHaplotypesInPopulation 6	
			echo "Sample" $sample "raw vcf generated and placed in" $rawvcf
		done < $sfile
		echo " Variant calling completed on all input files. Exiting.."

		echo ""
		echo "Step $answer completed on all input files in $samplelist. Now Exiting..."
		echo ""

elif [[ $answer == SnpFiltering ]] 
	then
		echo ""
		echo "you chose $answer, beginning..."
		echo ""
		#******# Perform GATK Operations
		while read sample; do
			echo "Variant call for" $rawvcf/$sample "found"
			# a. Breed3ck the variants that have been called by HaplotypeCaller using the VariantEval function. This will generate a table summarizing all types of variants
			java -jar $gatk -R $gallus -T VariantEval -eval $rawvcf/$sample".raw.vcf" --evalModule CountVariants -o $rawvcf/$sample".raw.VariantCounts" -noEV
			# b. Subset out only the SNP variants (using the SelectVariants tool) from the raw.vcf file. If you want to look at indels, a separate command can be used, where the selectType will be INDEL
			echo "Subsetting SNPs .. "
			java -jar $gatk -R $gallus -T SelectVariants --variant $rawvcf/$sample".raw.vcf" -o $snp/$sample".SNP.vcf" -selectType SNP
			# c. Let us apply some hard filtering to the new SNP file. The filters are applied using the --filterExpression. In this case, QD < 5 Fails the filter, and so on
			# Note that we are using the '||' or operator so it will test if one of them are failing. to use the and operater, it can be replaced with &&. 
			# The filter name has to be provided. Also the action will only annotate the VCF file. it will not remove them. that is a different step 
			java -jar $gatk -R $galus -T VariantFiltration -V $snp/$sample".SNP.vcf" --filterExpression "QD < 5.0 || FS > 55.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "HQFilter" -o $snp/$sample".HQSNP.vcf"
			# d. let us Breed3ck the contents of the HQSNP file. It is useful to compare the raw and the HQ files to see how the sample fared
			java -jar $gatk -R $gallus -T VariantEval -eval $snp/$sample".HQSNP.vcf" --evalModule VariantSummary -o $snp/$sample".HQSNP.VariantSummary" -noEV
			# e. subset a VCF file to select for PASS only variants, and retain header. 
			echo "Subsetting the passing SNPs only ..."
			awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $snp/$sample"HQSNP.vcf" > $snp/$sample".PASS.SNP.vcf"
			# f. Fun variant evaluation again to see how the counts compare in the passing set compared to raw
			java -jar $gatk -R $gallus -T VariantEval -eval $snp/$sapple".PASS.SNP.vcf" --evalModule CountVariants -o $snp/$sample".PASS.VariantCounts" -noEV
			echo "Variants successfully subset, and filtered"
		done < $sfile

		echo ""
		echo "Step $answer completed on all input files in $samplelist. Now Exiting..."
		echo ""

#******# SNP effect pipeline removed to preserve IP #******#

#******# VCFtools subsetting and comparison generation removed to preserve IP #******#

elif [[ $answer == Exit ]] 
	then
		echo ""
		echo "You chose the exit option. Good bye. Now exiting..."
		echo ""

else 
		echo ""
		echo "You have entered an unrecognized option. Please see logfile and/or README file to see correct entry format. Goodbye. "
		echo ""
fi





