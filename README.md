# cauris-data-flow 

cauris-data-flow is a set of python scripts for processing Candida auris sequencing data. This is meant to be used in combination with funQCD, nanofunQC, and nanofunsake.

To start, use process_new_samples.py and make_temp_reads_dir.py to add new samples to the metadata directory. This will also rename your reads for convenience.

Next, use setup.py and the temporary reads directory to run funQCD.

After running funQCD, the only step you should need to perform is running process_funqcd_shortread.py on the output directory:

```
python process_funqcd_shortread.py --input /scratch/esnitkin_root/esnitkin0/jjhale/funQCD/results/seqcoast_illumina_071526/
```

This will add your reads to the short_read and assembly datasets by default. The target directory can be changed with '-o'. 
Check the log files to ensure everything moved over correctly.





