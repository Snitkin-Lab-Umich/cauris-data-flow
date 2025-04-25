import pandas as pd
import subprocess
import os
import argparse

# the only move commands are currently commented out

def find_pass(final_qc_file):
    # take the final QC summary file and return the names of the samples that passed QC
    sample_pass,sample_fail = [],[]
    qc_df = pd.read_csv(final_qc_file, sep = '\t', header = 0)
    for i in range(qc_df.shape[0]):
        if qc_df.iloc[i]['QC_EVALUATION'] == 'PASS':
            sample_pass.append(qc_df.iloc[i]['Sample'])
        elif qc_df.iloc[i]['QC_EVALUATION'] == 'FAIL':
            sample_fail.append(qc_df.iloc[i]['Sample'])
    return([sample_pass,sample_fail])
    

def move_single(file1,source_dir,dest_dir,action='mv',debuglog = 'logs/debug.txt'):
    source_dir,dest_dir = [x+'/' if not x.endswith('/') else x for x in [source_dir,dest_dir]]
    with open(debuglog, 'a') as debug:
        if os.path.isfile(source_dir + file1) and os.path.isdir(dest_dir):
            subprocess.call([action,source_dir + file1,dest_dir + file1],stdout=debug,stderr=debug)
            _ = debug.write(' '.join([action,source_dir + file1,dest_dir + file1]) + '\n')
        else:
            _ = debug.write(f'ERROR: Failed to move {file1} to {dest_dir}\n')

def add_to_master(qc_file,master_file):
    with open(qc_file, 'r') as fhin, open(master_file, 'a') as fhout:
        next(fhin)
        for line in fhin:
            _ = fhout.write(line)
                
def init_master(qc_file,master_file):
    with open(qc_file, 'r') as fhin, open(master_file, 'w') as fhout:
        for line in fhin:
            _ = fhout.write(line)
            break
    print('Initialized master_qc_file.csv')


def move_all(batch_dir,qc_file,qcd_dir,debuglog = 'logs/debug.txt'):
    # take the directory of a specific batch and move the output files as needed
    # for samples that PASS:
    #   raw reads are moved to passed_samples under the batch directory and added to passed_samples.txt (no additional structure, and are seemingly moved out of raw_reads)
    #   trimmed reads go in clean_fastq_qc_pass_samples in the directory one level up (also no structure, these are moved out of the QCD trimmomatic directory)
    #   the QCD results we want to save are moved to similarly-named folders under assembly, with no batch names (i.e. spades/[sample_name]/[sample_name]_contigs_l1000.fasta and nothing else from spades)
    # for samples that FAIL:
    #   raw reads are moved to failed_samples under the batch directory and added to failed_samples.txt

    # change to directory, check paths, and make directories
    if not os.path.isdir(batch_dir):
        print(f'Could not locate batch directory at {batch_dir}')
        quit(1)
    batch_dir = os.path.abspath(batch_dir) + '/'
    os.chdir(batch_dir)
    if not os.path.isdir(qcd_dir):
        print(f'Could not locate QCD directory at {qcd_dir}')
        quit(1)
    qcd_dir = os.path.abspath(qcd_dir) + '/'
    for dname in ['passed_qc_samples','failed_qc_samples','logs','../clean_fastq_qc_pass_samples']:
        if not os.path.isdir(dname):
            subprocess.call(['mkdir',dname])
    with open(debuglog,'w') as fh:
        _ = fh.write(f'logs for {batch_dir}\n')

    # determine which samples and files passed qc
    sample_pass,sample_fail = find_pass(qc_file)

    # add this information to the master QC file
    if not os.path.isfile('../master_qc_summary.csv'):
        init_master(qc_file=qc_file,master_file='../master_qc_summary.csv')
    add_to_master(qc_file=qc_file,master_file='../master_qc_summary.csv')

    # move failed samples out of raw reads
    with open('failed_samples.txt','a') as failed_samples_log:
        for sample in sample_fail:
            filename = sample + '.fastq.gz'
            # move raw reads
            move_single(file1=filename,source_dir='raw_fastq/',dest_dir='failed_qc_samples/',debuglog=debuglog)
            _ = failed_samples_log.write(sample + '\n')
    
    # move passed samples out of raw reads
    # also move corresponding trimmed files to their own folder
    # also move qcd results to the assembly folder
    with open('passed_samples.txt','a') as passed_samples_log:
        for sample in sample_pass:
            # move raw reads
            filename = sample + '.fastq.gz'
            move_single(file1=filename,source_dir='raw_fastq/',dest_dir='passed_qc_samples/',debuglog=debuglog)
            _ = passed_samples_log.write(sample + '\n')

            # move trimmed reads, which are in qcd_dir/filtlong/[sample]/[sample].trimmed.fastq.gz
            filename_trim = sample + '.trimmed.fastq.gz'
            filt_dir = qcd_dir + 'filtlong/' + sample + '/'
            move_single(file1=filename_trim,source_dir=filt_dir,dest_dir='../clean_fastq_qc_pass_samples/',debuglog=debuglog)

            # move all qcd results
            # make a dict of all files to move
            # unfortunately, files and dirs have a variety of structures, so this is all manually entered
            actd = {}
            assembly_dir = '../../assembly/ONT/'

            # just the medaka assembly
            medaka_dir = qcd_dir + 'medaka/' + sample + '/' 
            medaka_file = sample + '_medaka.fasta'
            medaka_dest = assembly_dir + 'medaka/' + sample + '/'
            actd[medaka_file] = (medaka_dir,medaka_dest)

            # just the quast report
            quast_dir = qcd_dir + 'quast/' + sample + '/' + sample + '_medaka/'
            quast_file = 'report.txt'
            quast_dest = assembly_dir + 'quast_medaka/' + sample + '/'
            actd[quast_file] = (quast_dir,quast_dest)

            # the entire funannotate annotation results directory
            funannotate_dir = qcd_dir + 'funannotate/' + sample + '/annotate_results/'
            funannotate_dest = assembly_dir + 'funannotate/' + sample + '/annotate_results/'
            for file_name in os.listdir(funannotate_dir):
                funannotate_file = file_name
                actd[funannotate_file] = (funannotate_dir,funannotate_dest)

            busco_n_dir = qcd_dir + 'busco/busco_output_nucl/' + sample + '.scaffolds.fa/'
            busco_n_file = 'short_summary.specific.saccharomycetes_odb10.' + sample + '.scaffolds.fa.txt'
            busco_dest = assembly_dir + 'busco/' + sample + '/'
            actd[busco_n_file] = (busco_n_dir,busco_dest)

            busco_p_dir = qcd_dir + 'busco/busco_output_prot/' + sample + '.proteins.fa/'
            busco_p_file = 'short_summary.specific.saccharomycetes_odb10.' + sample + '.proteins.fa.txt'
            actd[busco_p_file] = (busco_p_dir,busco_dest)

            with open(debuglog,'a') as debug:
                for fname in actd:
                    source_dir = actd[fname][0]
                    dest_dir = actd[fname][1]
                    if not os.path.isdir(dest_dir):
                        subprocess.call(['mkdir','-p',dest_dir],stdout=debug,stderr=debug)
                        _ = debug.write(' '.join(['mkdir','-p',dest_dir]) + '\n')
                    move_single(file1=fname,source_dir=source_dir,dest_dir=dest_dir,action='cp',debuglog=debuglog)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--batch','-b',type=str,
        help='''Provide a path to the directory for this batch of Illumina sequencing data. 
        This should be in the project directory, under Sequence_data/ONT. Use an absolute path if possible.
        Example: /nfs/turbo/umms-esnitkin/Project_MDHHS_genomics/Sequence_data/ONT/2024-09-26_Plate1-to-Plate15/
        ''')
    parser.add_argument(
        '--qc_file','-q',type=str,
        help='''Provide a path to the final QC summary table for this batch of sequencing data. 
        This should be in the QCD directory, under [batch_name]/multiqc/[batch_name]_final_qc_summary.tsv. A relative path is OK here.
        ''')
    parser.add_argument(
        '--qcd_dir','-d',type=str,
        help='''Provide a path to the directory containing the QCD results. 
        This will usually just be the name of the batch. A relative path is OK here.
        ''')
    args = parser.parse_args()
    move_all(batch_dir=args.batch,qc_file=args.qc_file,qcd_dir=args.qcd_dir,debuglog = 'logs/debug.txt')

if __name__ == "__main__":
    main()







