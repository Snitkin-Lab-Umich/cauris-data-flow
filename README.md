# cauris-data-flow 

cauris-data-flow is a set of python scripts for processing Candida auris sequencing data. This is meant to be used in combination with funQCD, nanofunQC, and nanofunsake.

## Summary of steps

Each script is meant to be run in order. These primarily move files and create directories rather than performing any analysis. Each script has its own detailed help page, with a short summary given here.

### 1) rename_and_move.py

This is the first script to run, which uses a lookup table to move raw sequencing data to your project directory. It also gives each sample a simple, standardized name from the table and creates a new batch name.

Example usage:

```
python rename_and_move.py --input /nfs/turbo/umms-esnitkin/Project_Cauris/Sequence_data/metadata/sample_lookup/UM_illumina_sample_lookup.txt \
--project /nfs/turbo/umms-esnitkin/Project_Cauris/ --batch 2025_01_29_UM_illumina --type shortread --debug logs/UM_illumina_debug.txt
```

### 2) Run QCD

After your raw data files have been moved, run the appropriate QCD pipeline (such as funQCD). Point to the project directory as the location of the raw reads. Make sure the pipeline runs in its entirety without any errors. 

### 3) cleanup_move_qcd.py

Next, run this script. This moves the QCD results to the project directory and removes several very large intermediate files. Your QCD results should appear in your project directory next to the raw reads.

Example usage:

```
python cleanup_move_qcd.py -i /scratch/esnitkin_root/esnitkin0/jjhale/2025_01_29_pipline_test_v1/funQCD/results/2025-01-29_UM_illumina/ \
-o /nfs/turbo/umms-esnitkin/Project_Cauris/Sequence_data/illumina_fastq/2025_01_29_UM_illumina/ -n 2025_01_29_UM_illumina_QCD_Results
```

### 4) move_files_shortread.py (or move_files_longread.py)

This last script will examine the final QC results to determine which samples can be used for further analysis. Both the raw reads and the QCD outputs (such as protein annotations) will be sorted into different directories based on this information. Since there are different outputs and QC metrics for different data types, "move_files_shortread.py" should be used for Illumina data and "move_files_longread.py" should be used for Nanopore data.

Example usage:

```
python move_files_shortread.py -b /nfs/turbo/umms-esnitkin/Project_Cauris/Sequence_data/illumina_fastq/2025_01_29_UM_illumina/ \
-q 2025_01_29_UM_illumina_QCD_Results/multiqc/2025-01-29_UM_illumina_final_qc_summary.tsv \
-d 2025_01_29_UM_illumina_QCD_Results
```

After this script finishes, your files should have all been automatically moved to their appropriate directories.
