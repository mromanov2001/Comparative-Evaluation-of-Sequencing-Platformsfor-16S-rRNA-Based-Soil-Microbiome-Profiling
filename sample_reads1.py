import os
import random
from Bio import SeqIO
import gzip

random.seed(10)

def select_random_reads(input_file, output_file, num_reads=10000):
    # open .fastq.gz file
    with gzip.open(input_file, "rt") as handle:
        # Read all sequences from file
        reads = list(SeqIO.parse(handle, "fastq"))
    
    # if file has less than 10 000 reads, then show warning 
    if len(reads) < num_reads:
        print(f"Warning: File {input_file} contains less than {num_reads} reads. Selecting all reads.")
        num_reads = len(reads)
    
    # randomly choose num_reads reads
    selected_reads = random.sample(reads, num_reads)
    
    # save selected reads in a new .fastq file
    with gzip.open(output_file, "wt") as output_handle:
        SeqIO.write(selected_reads, output_handle, "fastq")

def process_directory(input_dir, output_dir, num_reads=10000):
    # check presence of the out dir, if it is absent tnen create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # process all files in dir
    for filename in os.listdir(input_dir):
        if filename.endswith(".fastq.gz"):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, f"{filename}")
            
            print(f"Processing {input_file}...")
            select_random_reads(input_file, output_file, num_reads=num_reads)
            print(f"Saved sampled reads to {output_file}")

# in и out dirs:
input_dir = "/path/to/input_dir"  #  dir with initial files

for reads in 10, 20, 25, 35:
  output_dir = f"/path/to/{reads}k_out_dir"  #  dir for new files with sampled reads
  process_directory(input_dir, output_dir, num_reads=reads*1000)
