import pandas as pd
import subprocess
import os
import argparse

def find_pass(final_qc_file):
    # take the final QC summary file and return the names of the samples that passed QC
    sample_pass,sample_fail = [],[]
    qc_df = pd.read_csv(final_qc_file, sep = '\t', header = 0)
    for i in range(qc_df.shape[0]):
        if qc_df.iloc[i]['QC_EVALUATION'] == 'PASS':
            sample_pass.append((qc_df.iloc[i]['Sample'],qc_df.iloc[i]['Filename']))
        elif qc_df.iloc[i]['QC_EVALUATION'] == 'FAIL':
            sample_fail.append((qc_df.iloc[i]['Sample'],qc_df.iloc[i]['Filename']))
    return([sample_pass,sample_fail])

def update_master_qc(master_qc_file, new_qc_file):
    # update a master QC summary file with the results from a new batch of samples
    # return the names of the samples that passed and failed QC in the new batch
    if not os.path.isfile(master_qc_file):
        print(f'Warning: {master_qc_file} does not exist. Making a new one with the contents of {new_qc_file}.')
        subprocess.call(['cp',new_qc_file,master_qc_file])
    else:
        # combine the two dataframes and overwrite the original master qc file
        # make a backup of the original master qc file as well
        subprocess.call(['cp',master_qc_file, master_qc_file.replace('.csv','_backup.csv')])
        master_df = pd.read_csv(master_qc_file, sep='\t')
        new_df = pd.read_csv(new_qc_file, sep='\t')
        master_df = pd.concat([master_df, new_df], ignore_index=True)
        master_df.to_csv(master_qc_file, sep='\t', index=False)
    sample_pass,sample_fail = find_pass(new_qc_file)
    return(sample_pass,sample_fail)

def add_funqcd_to_master(funqcd_dir, master_dir, master_qc_filename = 'master_qc_summary.csv'):
    # before doing anything, check if the columns in the new qc file match the names in the master qc file (if it exists)
    funqcd_batch_name = os.path.basename(funqcd_dir.rstrip('/'))
    new_qc_filename = None
    for f in os.listdir(os.path.join(funqcd_dir, 'multiqc')):
        if f.endswith('_final_qc_summary.tsv'):
            new_qc_filename = f
            break
    if new_qc_filename is None:
        print(f'Error: could not find the final QC summary file in {funqcd_dir}/multiqc')
        quit(1)

    master_qc_file = os.path.join(master_dir, master_qc_filename)
    if os.path.isfile(master_qc_file):
        master_qc_df = pd.read_csv(master_qc_file, sep='\t')
        funqcd_qc_df = pd.read_csv(os.path.join(funqcd_dir, 'multiqc', new_qc_filename), sep='\t')
        if not (len(master_qc_df.columns) == len(funqcd_qc_df.columns) and all(master_qc_df.columns == funqcd_qc_df.columns)):
            print(f'Error: The columns in {new_qc_filename} do not match those in {master_qc_filename}. Cannot add this batch to the master QC file.')
            print(f'Columns in {master_qc_filename}: {master_qc_df.columns.tolist()}')
            print(f'Columns in {new_qc_filename}: {funqcd_qc_df.columns.tolist()}')
            quit(1)
    
    # take a directory of funQCD outputs and update the dataset present in master_dir
    # move the funqcd output to master_dir
    funqcd_dest_dir = os.path.join(master_dir, funqcd_batch_name)
    funqcd_dest = os.path.join(funqcd_dest_dir, funqcd_batch_name + '_QCD_Results')
    if os.path.isdir(funqcd_dest_dir):
        print(f'Error: A directory named {funqcd_batch_name} already exists in {master_dir}')
        quit(1)
    subprocess.call(['mkdir','-p',os.path.join(funqcd_dest_dir,'logs')])
    debuglog = os.path.join(funqcd_dest_dir,'logs','debuglog.txt')
    with open(debuglog,'a') as debug:
        _ = debug.write(f'logs for {funqcd_batch_name}\n')
        # exclude specific intermediate files and directories when moving funqcd_dir to funqcd_dest
        exclude_list = ['run_saccharomycetes_odb10','K21','K33','K55','K77','pipeline_state','tmp','annotate_misc','repeatmasker']
        cmd1 = ['rsync','-r']
        for exclude in exclude_list:
            cmd1.append('--exclude=' + exclude)
        cmd1+=[funqcd_dir,funqcd_dest]
        subprocess.call(cmd1)
        _ = debug.write(' '.join(cmd1) + '\n')
    # identify the location of the QC summary file in funqcd_dest
    qc_dir = os.path.join(funqcd_dest,'multiqc')
    qc_file = None
    for f in os.listdir(qc_dir):
        if f.endswith('_final_qc_summary.tsv'):
            qc_file = os.path.join(qc_dir,f)
            break
    if qc_file is None:
        print(f'Error: could not find the final QC summary file in {qc_dir}')
        quit(1)
    new_qc_file = os.path.join(funqcd_dest,'multiqc',new_qc_filename)
    # use this to update the master QC summary file, as well as identifying which samples passed and failed QC
    passed_samples, failed_samples = update_master_qc(master_qc_file=master_qc_file, new_qc_file=new_qc_file)

    # record all failed samples in a text file
    failed_samples_file = os.path.join(funqcd_dest_dir,'failed_samples.txt')
    with open(failed_samples_file,'w') as fh:
        for sample,filename in failed_samples:
            _ = fh.write(sample + '\n')
    
    # record all passed samples in a text file
    # in addition, move the corresponding trimmed reads to pass_trimmed_reads in master_dir
    passed_samples_file = os.path.join(funqcd_dest_dir,'passed_samples.txt')
    with open(passed_samples_file,'w') as fh:
        for sample, filename in passed_samples:
            _ = fh.write(sample + '\n')
            # move trimmed reads, which are in trimmomatic/[sample]/[sample]_R*_trim_paired.fastq.gz
            filename_trim_r1,filename_trim_r2 = sample + '_R1_trim_paired.fastq.gz', sample + '_R2_trim_paired.fastq.gz'
            trim_r1 = os.path.join(funqcd_dest,'trimmomatic',sample,filename_trim_r1)
            trim_r2 = os.path.join(funqcd_dest,'trimmomatic',sample,filename_trim_r2)
            if not os.path.isfile(trim_r1) or not os.path.isfile(trim_r2):
                print(f'WARNING: Could not locate or move trimmed reads for {sample}!! Expected files: {trim_r1} and {trim_r2}')
            else:
                cmd1 = ['mv',trim_r1,os.path.join(master_dir,'pass_trimmed_reads')]
                cmd2 = ['mv',trim_r2,os.path.join(master_dir,'pass_trimmed_reads')]
                subprocess.call(cmd1)
                subprocess.call(cmd2)

            # move all qcd results
            # make a dict of all files to move
            # unfortunately, files and dirs have a variety of structures, so this is all manually entered
            actd = {}
            assembly_dir = os.path.join(os.path.dirname(master_dir),'assembly','illumina')
            if not os.path.isdir(assembly_dir):
                print(f'Error: could not locate assembly directory at {assembly_dir}')
                quit(1)

            # just the spades assembly
            spades_src = os.path.join(funqcd_dest,'spades',sample,f'{sample}_contigs_l1000.fasta')
            if not os.path.isfile(spades_src):
                print(f'Warning: could not locate spades assembly for {sample} at {spades_src}')
            else:
                spades_dest = os.path.join(assembly_dir,'spades',sample,f'{sample}_contigs_l1000.fasta')
                actd['spades'] = (spades_src,spades_dest)

            # just the quast report
            quast_src = os.path.join(funqcd_dest,'quast',sample,'report.txt')
            if not os.path.isfile(quast_src):
                print(f'Warning: could not locate quast report for {sample} at {quast_src}')
            else:
                quast_dest = os.path.join(assembly_dir,'quast',sample,'report.txt')
                actd['quast'] = (quast_src,quast_dest)

            # the entire funannotate annotation results directory
            funannotate_src = os.path.join(funqcd_dest,'funannotate',sample,'annotate_results')
            if not os.path.isdir(funannotate_src):
                print(f'Warning: could not locate funannotate results for {sample} at {funannotate_src}')
            else:
                funannotate_dest = os.path.join(assembly_dir,'funannotate',sample,'annotate_results')
                actd['funannotate'] = (funannotate_src,funannotate_dest)

            # the busco nucleotide results
            busco_n_file = 'short_summary.specific.saccharomycetes_odb10.' + sample + '.scaffolds.fa.txt'
            busco_n_src = os.path.join(funqcd_dest,'busco','busco_output_nucl',f'{sample}.scaffolds.fa',busco_n_file)
            if not os.path.isfile(busco_n_src):
                print(f'Warning: could not locate busco nucleotide results for {sample} at {busco_n_src}')
            else:
                busco_n_dest = os.path.join(assembly_dir,'busco',sample,busco_n_file)
                actd['busco_n'] = (busco_n_src,busco_n_dest)

            # the busco protein results
            busco_p_file = 'short_summary.specific.saccharomycetes_odb10.' + sample + '.proteins.fa.txt'
            busco_p_src = os.path.join(funqcd_dest,'busco','busco_output_prot',f'{sample}.proteins.fa',busco_p_file)
            if not os.path.isfile(busco_p_src):
                print(f'Warning: could not locate busco protein results for {sample} at {busco_p_src}')
            else:
                busco_p_dest = os.path.join(assembly_dir,'busco',sample,busco_p_file)
                actd['busco_p'] = (busco_p_src,busco_p_dest)

            with open(debuglog,'a') as debug:
                for fname in actd:
                    src = actd[fname][0]
                    dest = actd[fname][1]
                    if not os.path.isdir(os.path.dirname(dest)):
                        subprocess.call(['mkdir','-p',os.path.dirname(dest)])
                        debug.write(f'made directory {os.path.dirname(dest)}\n')
                    cmd = ['cp','-r',src,dest]
                    debug.write(' '.join(cmd) + '\n')
                    subprocess.call(cmd,stdout=debug,stderr=debug)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to the input directory. This should contain the funQCD results.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a path output directory. This should be the dataset you want to update with new results. The master QC summary file should be located here.''',
        default='/nfs/turbo/umms-esnitkin/Project_Cauris/Sequence_data/short_read'
        )
    args = parser.parse_args()
    add_funqcd_to_master(funqcd_dir=args.input, master_dir=args.output)

if __name__ == "__main__":
    main()







