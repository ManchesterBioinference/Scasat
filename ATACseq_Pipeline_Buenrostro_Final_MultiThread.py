
# coding: utf-8

# # ATACseq pipeline

# This jupyter-notebook automates the ATACseq pre-processing pipeline as well as the basic analysis

# In the following section we assume that:
# - The reference genome has to be already indexed with bowtie
# - The tutorials are introduced with the example file mentioned in the tutorial

# ## Importing python modules


import subprocess, os, sys, signal
import rpy2
from joblib import Parallel, delayed
import multiprocessing
import threading
import time
from timeit import default_timer as timer


# ## Setting the folders and tools

# First we configure the __Folders__


# Setting the output folder where all the results would be stored
outputFolder = '/../../mqbsxsm2/scratch/Buenrostro_pipeline_Final/Human_data_analysis/output/'

# inputFolder : Folder name with all the fastq files
intputFolder = '/../single-cell/users/mqbsxsm2/Buenrostro_ATAC-seq_data/Homosapiens_fastq_files_selected'

# Blacklisted sequences
blackListFile = '/../../mqbsxsm2/Packages/Blacklist_regions/consensusBlacklist.bed'

# Setting up the genome folder
ref_genome = '/../../mqbsxsm2/Bowtie2_Human_genome/genome'


# ### Setting up the softwares
# 
# We configure the __software__ parameters here. If the required softwares are not in the PATH you can manually set their location here

bowtie_path = 'module load apps/bowtie2/2.2.6/gcc-4.8.5; bowtie2'
samtools_path = 'module load apps/samtools/1.2/gcc-4.8.5; samtools'
macs2_path = 'module load apps/macs/2.1.1.20160309/gcc-4.8.5+numpy-1.9.2; macs2'
bdg2bw_path = '/mnt/mr01-home01/mqbsxsm2/Packages/bdg2bw'
picard_path = 'module load apps/picard/2.1.0/bin; picard.jar'
intersectBed_path = 'module load apps/bedtools/2.25.0/gcc-4.8.5; intersectBed'
MarkDuplicate_path = 'module load apps/picard/2.1.0/bin; MarkDuplicates'
Trimmomatic_path = 'module load apps/trimmomatic/0.36/noarch; trimmomatic'
fastqc_path ='module load apps/fastqc/0.11.3/noarch; fastqc'
Trimmomatic_Adapter = '/opt/gridware/local/el7/pkg/apps/trimmomatic/0.36/noarch/share/adapters/TruSeq3-PE-2.fa:2:30:10'


# Defining the software variables that would be called later on for execution

softwares = {    
    'bowtie2': bowtie_path,
    'macs2': macs2_path,
    'bdg2bw': bdg2bw_path,
    'samtools': samtools_path,
    'picard': picard_path,
    'intersectBed': intersectBed_path,
    'MarkDuplicate': MarkDuplicate_path,
    'trimmomatic': Trimmomatic_path,
    'fastqc': fastqc_path,
    'Trimmomatic_Adapter': Trimmomatic_Adapter
}


# ### Creating the output folders
# We now create the output folders where all the processed files will be stored

output_folders = [ 'Bowtie_files', 'Blacklist_removed' # Trimmo, Bowtie, Blacklist Removed Files
                 , 'Trimmomatic_Files', 'Fastqc_SE_Files' 
                 , 'Duplicates_removed', 'Macs2_files'              # Synced R1 and R2
                 , 'Merged_BAM', 'Merged_Macs2_files'
                 , 'Bam_Files_Filtered_on_Library_Size', 'Merged_Filtered_BAM'
                 ,  'Merged_Filtered_Macs2_files', 'Filtered_Macs2_files'
                 ]


for folder in output_folders:
    if not os.path.exists(os.path.join(outputFolder, folder)):
        os.makedirs(os.path.join(outputFolder, folder))


# ### Custom functions

# Here we define __custom__ functions that we will use for file processing

remove_extension = lambda x: x.split('.')[0]


# The following functions deals with the inputs and the outputs


def get_args(num_threads, read1, read2, ref_genome, blackListFile, output_folders):
    '''Set the input and output path for a given pair of reads'''
    r1_shortname = remove_extension(os.path.basename(read1))
    sample_name = r1_shortname.split('_')[0] + '_'+ r1_shortname.split('_')[1]

    args = {
        'num_threads': num_threads,
        'r1_input': read1,
        'r2_input': read2,
        'ref_genome': ref_genome,
        'blackListFile': blackListFile,
        'sample_name': sample_name,
    }
    
    ## output_paths = {folder: os.path.join(outputFolder, folder, r1_shortname) for folder in output_folders}
    output_paths = {folder: os.path.join(outputFolder, folder, sample_name) for folder in output_folders}
    
    return dict(args, **output_paths)



cmds = [
    
    '{trimmomatic} PE -phred33 {r1_input} {r2_input} {Trimmomatic_Files}_r1_paired.fq {Trimmomatic_Files}_r1_unpaired.fq {Trimmomatic_Files}_r2_paired.fq {Trimmomatic_Files}_r2_upaired.fq ILLUMINACLIP:{Trimmomatic_Adapter}',
    
    '{bowtie2} -p 12 -X 2000 --dovetail -x {ref_genome} -1 {Trimmomatic_Files}_r1_paired.fq -2 {Trimmomatic_Files}_r2_paired.fq -S {Bowtie_files}.sam 2> {Bowtie_files}_allignment.log',
    
    #'{bowtie2} -p 12 -X 2000 --dovetail -x {ref_genome} -1 {r1_input} -2 {r2_input} -S {Bowtie_files}.sam 2> {Bowtie_files}_allignment.log',
    
    '{samtools} view -SbhF 4 -f2 -q30 {Bowtie_files}.sam > {Bowtie_files}.bam' ,
    
    '{intersectBed} -v -abam {Bowtie_files}.bam -b {blackListFile} > {Blacklist_removed}_noexclu.bam ' ,
    
    '{samtools} sort {Blacklist_removed}_noexclu.bam {Blacklist_removed}_noexclu_sorted' ,
    
    '{MarkDuplicate}  INPUT={Blacklist_removed}_noexclu_sorted.bam OUTPUT={Duplicates_removed}_nodup.bam METRICS_FILE={Duplicates_removed}_nodup_stats REMOVE_DUPLICATES=True' ,
    
    '{samtools} sort {Duplicates_removed}_nodup.bam {Duplicates_removed}_nodup_sorted' ,
    
    '{samtools} index {Duplicates_removed}_nodup_sorted.bam' ,
    
    '{samtools} view -b {Duplicates_removed}_nodup_sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX  > {Duplicates_removed}_nodup_sorted_cleaned.bam' ,
        
    '{samtools} index {Duplicates_removed}_nodup_sorted_cleaned.bam',
    
    '{samtools} idxstats {Duplicates_removed}_nodup_sorted_cleaned.bam > {Duplicates_removed}_nodup_sorted_cleaned_chrom_stat.txt' ,
    
    '{macs2} callpeak -t {Duplicates_removed}_nodup_sorted_cleaned.bam -n {Macs2_files} -p 0.0001 -g hs -f BAMPE --nomodel --nolambda -B --keep-dup all --call-summits',
        
    '{bdg2bw} {Macs2_files}_treat_pileup.bdg /mnt/mr01-home01/mqbsxsm2/Packages/Homer_4.6/bin/hg19.chrom.sizes',
        
]


# Get the reads from the `inputFolder`
for root, folders, files in os.walk(intputFolder):
    files = [f for f in files if not f.startswith('.')] #remove hidden files if there exist
    reads1 = sorted([os.path.join(root, f) for f in files if '_1' in f])
    reads2 = sorted([os.path.join(root, f) for f in files if '_2' in f])


# Print the read names. This is to make sure that the correct files are picked up. This can be run just once as it is expected to produce large output. 
##print (reads1, reads2)


# The platform is now set to run the commands. Below we run the command sequentially. However, it is always good to first check that the right commands are executed. To test that one can uncomment the print command


def run_cmd(cmds, args, sema):
    try:
        for cmd in cmds: 
            print threading.currentThread().getName(), "running the command: ", cmd.format(**args)
            subprocess.call(cmd.format(**args), shell=True)
    finally:
        # Thread has finished its work
        print threading.currentThread().getName(), "finished."
        # Inform the main thread that it can start a new thread
        sema.release()


# Looking at the number of cores in the node
## print("Total Cores:",multiprocessing.cpu_count())


# __Multithread starts __

# We restrict the number of threads python will start (each will num trimmomatic then STAR, ...)
# Each software app in those threads will start their own (java) threads. So restrict the number.
# num_python_threads * num_software_threads must equal $NSLOTS (the number of cores reserved in batch)
# Number of theads python should start
num_python_threads=4


# Number of threads each software (trimmomatic, STAR, ...) should start
num_software_threads=4


# Keep a list of thread objects we create - one per input file.
# But we only allow a small number of python threads to run at any one time.
thread_list = []
thread_pool_sema = threading.BoundedSemaphore(value=num_python_threads)


##print("Starting ", num_python_threads, "threads...", flush=True)
start_t = timer()
for read1, read2 in zip(reads1, reads2):
    args = get_args(num_software_threads, read1, read2, ref_genome, blackListFile,output_folders)
    args = dict(args, **softwares)
    t = threading.Thread(target=run_cmd, args = (cmds,args,thread_pool_sema))
    t.daemon = True       
    thread_list.append(t)
    # This will block if limit of running python threads reached.
    # Only when another thread finishes will we be release to continue.
    thread_pool_sema.acquire()
    try:
        t.start()
        print "Main thread has started", t.name
    except:
        print "Main thread cannot start thread", t.name

# Wait here until all items have been processed
for t in thread_list:
    t.join(timeout=None)
print "Threaded tasks done in ", timer()-start_t, "seconds"



