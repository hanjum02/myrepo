# Overview for Track 1 Pipeline
   The purpose for this code is to perform differential expression between samples starting from a desired SRA file. The SRA file is then converted to paired-end fastq files which will then have their TPM quanitified using kallisto and BioPython. The genome's number of coding sequences, as well as the min, median, mean, and max TPM will be calculated and written to the log file. The output from kallisto is then used to pass through R's sleuth package to find differentially expressed genes between two
   conditions. The most differentially expressed CDS/protein fasta file is retrieved and used as a blast+ input to search the nucleotide database limited to the Betaherpesvirinae subfamily. This will return the subject accession, percent identity, alignment length, start and end of alignment in query, start and end of alignment in subject, bitscore, e-value, and subject title. 

# Dependencies
 - os (used to run command line commands from Python)
 - sratoolkit (used for converting SRA files to paired-end fastq files)
 - kallisto (used for building transcriptome index for HCMV and quantifying the TPM statistics for each CDS in each transcriptome)
 - statistics (used for calculating statistical values such as mean, median, min, and max)
 - sleuth (used for differential expression analysis in R)
 - dplyr (used for differential expression analysis in R)
 - blast+ (used to run NCBI blast)
 - BioPython (used for file parsing, generating input for kallisto, and biological calculations)  
 - pathlib (used for finding, establishing, and accessing files and directories from the user's system)
 - argparse (used for parsing command line functions, arguments, and subarguments)

# Step 1 Process
  1) Make sure the SRA-Toolkit package (sratoolkit) is installed since this will be used to convert the SRA file to paired-end fastq files
  2) Start by collecting the SRA accession/SRR numbers (in the format of SRR____) for each sample by visiting the respective NCBI page and looking under the "Runs" subheader
  3) Use the prefetch SRA-Toolkit function followed by the SRA accession/SRR number, this downloads the SRA file to your system. Example below:
      example_user:~$ prefetch SRR5660030
  4) Now that the SRA file is downloaded, use the fasterq-dump SRA-Toolkit function followed by the SRA accession/SRR number to convert it to a paired end fastq file. Use
     the --split-files & --skip-technicals options to create separate files for each paired-end read and skip technical reads respectively. After executing this function,
     there should now be two fastq files in your directory where they can be moved or analyzed as needed. Example below:
      example_user:~$ fasterq-dump SRR5660030 --split-files --skip-technical

# How to Run the Code
- To run the code, ensure that each necessary file is downloaded to your system (like the SRRs obtained from Step 1 as well as the sleuth.r R script required to run
  differential expression analysis), and more specifically, to your current directory where you are running the code from. Then, be sure to run the code by using the command line to 
  call python followed by the path to where the PipelineProject script is located and then finally including the name of the directory where you would like the results stored. Example below:
        example_user::~$ python "PathtoPipelineCode/PipelineCode.py" Output_Directory