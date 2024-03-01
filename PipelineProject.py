# imports all necessary packages to perform specified functions 
from Bio import Entrez, SeqIO
import os
import statistics
from pathlib import Path
import argparse 

def parse_arguments():
    parser = argparse.ArgumentParser(description='Track 1 Differential Expression Pipeline')
    parser.add_argument('output_dir', help='Output directory path', type=str)
    return parser.parse_args()

# defines function that acquires the HCMV genome using HCMV accesion number and Biopython's Entrez & SeqIO functions
def get_HCMVgenome(accessionID):
    Entrez.email = "hanjum@luc.edu"
    # searches nucleotide database using given HCMV accessionID in the GenBank file format
    handle = Entrez.efetch(db = "nucleotide", id = accessionID, rettype = "gb", retmode = "text")
    # variable stores the result of SeqIO function's reading of the GenBank file 
    record = SeqIO.read(handle, "genbank")
    return record

# defines function to extract all CDS found in the GenBank file/record retrieved from the previous function and write them to a file
def get_CDS(record, outfile):
    # opens output file that will serve as a fasta file to write the CDS features from the GenBank file/record to 
    with open(outfile, "w") as file:
        # cumulative sum variable that tallies the number of CDS found in the record
        CDS = 0
        # iterates through all features of the GenBank record to find all CDS and if the current feature is a CDS, increment CDS by 1
        for feature in record.features:
            if feature.type == "CDS":
                CDS += 1
                # extracts the corresponding protein ID for the current CDS using .qualifiers SeqFeature function
                proteinID = feature.qualifiers["protein_id"][0]
                # variable stores the sequence for the current CDS using .extract and .seq
                sequence = feature.extract(record.seq)
                # writes all the CDS sequences and their protein IDs to output file in the format of a fasta file
                file.write(">" + proteinID + "\n" + str(sequence) + "\n")
    return CDS

# defines function that runs the kallisto indexing command using os.system to run it on the command line from Python
def transcriptome_index(fasta, index):
    index_str = str(index)
    fasta_str = str(fasta)
    # sets the string concatenation for the kallisto index command that will be passed through the os.system() function 
    Indexcommand = "kallisto index -i " + index_str + " " + fasta_str
    os.system(Indexcommand)

# defines function that runs kallisto quantification command using os.system to run it on the command line from Python
def TPM_quantification(index, sampleID, fastq1, fastq2, output_dir):
    # sets the string which will serve as the output directory for each sample, storing the results of each kallisto quant command 
    output_directory = output_dir / f"TPM_Quantification_{sampleID}"
    output_directory.mkdir(parents=True, exist_ok=True)
    # sets the string for the kallisto quant command that will be passed through the os.system() function
    TPMcommand = f"kallisto quant -i {index} -o {output_dir} -b 30 -t 4 {fastq1} {fastq2}"
    os.system(TPMcommand)
    # Move the abundance.tsv file to the correct subdirectory
    abundance_file = output_dir / "abundance.tsv"
    new_abundance_file = output_directory / "abundance.tsv"
    os.rename(abundance_file, new_abundance_file)

# defines function which calculates the minimum, median, mean, and maximum TPM from the results in the abundance.tsv
def TPM_statistics(output_dir):
    # using os.path.join function, the abundance.tsv file can be pointed to and accessed
    abundance = output_dir / "abundance.tsv"
    TPM_values = []
    # opening and reading the abundance.tsv file, the TPM results are extracted 
    with open(abundance, "r") as file:
        next(file)
        for line in file:
            # splits the file's contents based on "\t", aka tab-delimiter and stores in list format
            fields = line.strip().split("\t")
            # appends the TPM value stored in the corresponding element to the TPM_values list
            TPM_values.append(float(fields[4]))
    # calculates minimum TPM value for the sample
    TPM_min = min(TPM_values)
    # calculates median TPM value for the sample
    TPM_median = statistics.median(TPM_values)
    # calculates mean TPM value for the sample
    TPM_mean = statistics.mean(TPM_values)
    # calculates maximum TPM value for the sample
    TPM_max = max(TPM_values)
    # returns each statistical value and rounds their values to the thousandths 
    return str(round(TPM_min, 3)), str(round(TPM_median, 3)), str(round(TPM_mean, 3)), str(round(TPM_max, 3))

# function that looks at the SRR samples then creates directories for each sample to hold their respective abundance.tsv files
def process_SRRs(samples, output_dir, index, fastqs):
    for sampleID, sample_info in samples.items():
        condition = sample_info["condition"]
        fastq1, fastq2 = sample_info["fastq_pair"]
        TPM_directory = output_dir / f"TPM_Quantification_{sampleID}"
        TPM_directory.mkdir(parents = True, exist_ok = True)
        TPM_quantification(index, sampleID, fastq1, fastq2, output_dir)
        min_tpm, med_tpm, mean_tpm, max_tpm = TPM_statistics(TPM_directory)

# function that retrieves the most differentially expressed CDS 
def get_mostDE_CDS(results_file):
    with open(results_file, 'r') as file:
        next(file)
        line = file.readline().strip()
        DE_CDS_id = line.split()[0]
    return DE_CDS_id

# function that retrieves the corresponding fasta file for the most DE CDS using BioPython
def get_mostDE_fasta(CDS_id):
    Entrez.email = "hanjum@luc.edu"
    # searches nucleotide database using given HCMV accessionID in the GenBank file format
    handle = Entrez.efetch(db = "protein", id = CDS_id, rettype = "fasta", retmode = "text")
    # variable stores the result of SeqIO function's reading of the GenBank file 
    record = SeqIO.read(handle, "fasta")
    return record

# function that downloads the NCBI sequences using 'datasets downlaod' 
def download_HCMVseqs(searchterm):
    search_command = "datasets download virus genome taxon " + searchterm + " --refseq --include genome"
    os.system(search_command)
    unzipcommand = "unzip ncbi_dataset.zip"
    os.system(unzipcommand)

# function that creates a local databse using a given databaseID and a fasta file to use as input 
def create_localdb(file_input, database):
    # Convert PosixPath objects to strings
    file_input_str = str(file_input)
    database_str = str(database)
    makeblast_command = "makeblastdb -in " + file_input_str + " -out " + database_str + " -title " + database_str + " -dbtype nucl"
    os.system(makeblast_command)

# function that runs the blast+ command from the command line using os.system and a string which holds the command to use
def run_blast(query_fasta, database, out_put):
    blast_command = "tblastn -query " + query_fasta + " -db " + database + " -out " + out_put + " -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' -max_target_seqs 10"
    os.system(blast_command)

def main():
    
    # using pathlib and argparse to determine where to set up necessary files and directories for access later on
    args = parse_arguments()
    output_dir = Path(args.output_dir)
    Path(output_dir).mkdir(parents = True, exist_ok = True)
    nuc_seqs = Path("ncbi_dataset/data/genomic.fna")
    output_file = output_dir / "PipelineProject.log"
    r_script = Path("mycode.R")
    r_results = Path("fdr05_results.txt")

    fastqs = list(Path("PipelineFastqs").glob("*.fastq"))

    hcmv_acc = "NC_006273.2"
    hcmv_genome = get_HCMVgenome(hcmv_acc)

    CDS_file = output_dir / "CDS.fasta"
    CDS_count = get_CDS(hcmv_genome, CDS_file)

    indexID = output_dir / "HCMV_index.idx"
    transcriptome_index(CDS_file, indexID)

    # dictionary holding the conditions and fastqs of the SRR samples
    samples = {
        "SRR5660030": {
            "condition": "2dpi",
            "fastq_pair": ("PipelineFastqs/SRR5660030_1.fastq", "PipelineFastqs/SRR5660030_2.fastq")},
        "SRR5660033": {
            "condition": "6dpi",
            "fastq_pair": ("PipelineFastqs/SRR5660033_1.fastq", "PipelineFastqs/SRR5660033_2.fastq")},
        "SRR5660044": {
            "condition": "2dpi",
            "fastq_pair": ("PipelineFastqs/SRR5660044_1.fastq", "PipelineFastqs/SRR5660044_2.fastq")},
        "SRR5660045": {
            "condition": "6dpi",
            "fastq_pair": ("PipelineFastqs/SRR5660045_1.fastq", "PipelineFastqs/SRR5660045_2.fastq")}
    }

    # Glob for FASTQ files in the directory and sort them for paired-end reads
    fastqs = sorted(Path("PipelineFastqs").glob("*.fastq"))

    # Group paired-end reads together
    fastq_pairs = [(fastqs[i], fastqs[i+1]) for i in range(0, len(fastqs), 2)]

    process_SRRs(samples, output_dir, indexID, fastqs)

    # opens and writes to the log file all the necessary information gathered from running TPM quant and statistics functions
    with open(output_file, "a") as file:
        file.write("The HCMV genome " + hcmv_acc + " has " + str(CDS_count) + " CDS.\n")
        file.write("{:<12}\t{:>8}\t{:>8}\t{:>8}\t{:>8}\t{:>8}\n".format("sample", "condition", "min_tpm", "med_tpm", "mean_tpm", "max_tpm"))
        for sampleID, condition_info in samples.items():
            condition = condition_info["condition"]
            TPM_directory = output_dir / f"TPM_Quantification_{sampleID}"
            min_tpm, med_tpm, mean_tpm, max_tpm = TPM_statistics(TPM_directory)
            file.write("{:<12}\t{:>8}\t{:>8}\t{:>8}\t{:>8}\t{:>8}\n".format(sampleID, condition, min_tpm, med_tpm, mean_tpm, max_tpm))

    # runs the R script from the command line
    os.system("Rscript " + str(r_script))

    # opens and writes the headers for the sleuth r script DE analysis values
    with open(output_file, "a") as file:
        file.write("{:<12}\t{:<20}\t{:<20}\t{:<20}\n".format("target_id", "test_stat", "pval", "qval"))
        
    # opens and reads the results of running the R script, then writes the results to the log file
    with open(r_results, "r") as results_file:
        next(results_file)  # Skip header line
        for line in results_file:
            info = line.strip().split()
            target_id = info[0]
            pval = info[1]
            qval = info[2]
            test_stat = info[3]
            with open(output_file, "a") as log_file:
                log_file.write("{:<12}\t{:<20}\t{:<20}\t{:<20}\n".format(target_id, test_stat, pval, qval))

    mostDE_CDS = get_mostDE_CDS(r_results)

    DE_fasta = get_mostDE_fasta(mostDE_CDS)

    # sets up the necessary files needed to run the download_HCMVseqs() and createlocaldb() functions
    searchterm = "Betaherpesvirinae"
    download_HCMVseqs(searchterm)
    create_localdb(nuc_seqs, searchterm)

    # sets up the necessary files needed to run the blast+ search as called by the run_blast() function
    SeqIO.write(DE_fasta, "query.fasta", "fasta")
    query = "query.fasta"
    blast_results = "HCMV_tblastn_results.csv"
    run_blast(query, searchterm, blast_results)

    # opens and reads the results of the blast+ search, then writes the top 10 hits to the log file
    with open(blast_results, 'r') as file:
        with open(output_file, 'a') as log:
            log.write("{:<12}\t{:>8}\t{:>6}\t{:>6}\t{:>6}\t{:>8}\t{:>8}\t{:>8}\t{:>10}\t{}\n".format("sacc", "pident", "length", "qstart", "qend", "sstart", "send", "bitscore", "evalue", "stitle"))
            # establishes a cumulative sum variable to track the number of lines processed
            line_count = 0  
            for line in file:
                # checks for only the top 10 hits 
                if line_count >= 10:  
                    # once the top 10 hits have been checked, exits the loop
                    break  
                # strips and splits the necessary values that will be printed
                values = line.strip().split('\t')
                formatted_line = "{:<12}\t{:>8}\t{:>6}\t{:>6}\t{:>6}\t{:>8}\t{:>8}\t{:>8}\t{:>10}\t{}\n".format(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9])
                log.write(formatted_line)
                # cumulative sum variable increments by 1 for each line iterated
                line_count += 1

if __name__ == "__main__":
    main()