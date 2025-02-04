import pandas as pd
import os
import subprocess
import glob
import argparse

# this file needs a csv with a standard format:
# [sample_id] [new_name] [path_to_reads]

def move_and_rename(input_path,output_path,data_type,debuglog='logs/debuglog.txt'):
    # this function takes a three-column csv and uses it to move data to a new directory with a new name
    suffixes = ('.fastq.gz','.fq.gz')
    # read conversion table
    conv = pd.read_csv(input_path, sep = '\t', header = 0)
    for i in range(conv.shape[0]):
        if data_type == 'longread':
            original_file = conv.iloc[i][2]
            if not original_file.endswith(suffixes):
                print(f'Unexpected file format in {original_file}')
                quit(1)
            new_file = output_path + conv.iloc[i][1] + '.fastq.gz'
            with open(debuglog,'a') as debug:
                command = ['cp',original_file,new_file]
                subprocess.call(command,stdout=debug,stderr=debug)
                _ = debug.write(' '.join(command) + '\n')
        if data_type == 'shortread':
            search_str = conv.iloc[i][2].replace('*','_*')
            # the conversion files have file names where a * immediately follows a number, which creates ambiguity between 1 and 10 and such
            # ex: ONT_Sequencing/Lojek_D5A_polished_results/Lojek_D5A_1*.fastq.gz and ONT_Sequencing/Lojek_D5A_polished_results/Lojek_D5A_10*.fastq.gz
            # this search_str line will need to be changed if the conversion files are formatted differently
            flist = sorted(glob.glob(search_str))
            if len(flist) != 2:
                print(f'Error while searching for files at path {conv.iloc[i][2]}')
                quit(1)
            original_file_r1,original_file_r2 = sorted(glob.glob(search_str))
            # file paths for illumina data will contain * where R1 and R2 would be, in addition to other idiosyncracies
            # sorting the list should ensure that R1 and R2 are assigned correctly
            if '_R1' not in original_file_r1 or '_R2' not in original_file_r2 or not (original_file_r1.endswith(suffixes) and original_file_r2.endswith(suffixes)):
                print(f'Error assigning R1 and R2 for file path {conv.iloc[i][2]}')
                quit(1)
            new_file_r1 = output_path + conv.iloc[i][1] + '_R1.fastq.gz'
            new_file_r2 = output_path + conv.iloc[i][1] + '_R2.fastq.gz'
            with open(debuglog,'a') as debug:
                command1,command2 = ['cp',original_file_r1,new_file_r1],['cp',original_file_r2,new_file_r2]
                subprocess.call(command1,stdout=debug,stderr=debug)
                subprocess.call(command2,stdout=debug,stderr=debug)
                _ = debug.write(' '.join(command1) + '\n')
                _ = debug.write(' '.join(command2) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to a conversion file. This should be a three-column CSV consisting of SAMPLE_ID, NEW_ID, and PATH_TO_READS (in that order).
        The file(s) at PATH_TO_READS will be given a new name based on NEW_ID and moved to the directory specified by --project and --batch.
        '''
        )
    parser.add_argument(
        '--project','-p',type=str,
        help='''Provide a path to the project directory. This is the directory containing Sequence_data and Analysis.'''
        )
    parser.add_argument(
        '--batch','-b',type=str,
        help='''Provide a name for this batch of samples. A separate directory will be created with this name. Prefix this name with the date if possible.'''
        )
    parser.add_argument(
        '--type','-t',type=str,choices=['shortread','s','illumina','longread','l','nanopore','ont'],
        help='''Provide the type of sequencing data (short-read or long-read).'''
        )
    parser.add_argument(
        '--qc','-q',type=str,
        help='''(Optional) Provide a path to one of the QC and annotation pipelines (such as funQCD or nanofunQC). 
        This should be the path to the main directory containing results, workflow, config, etc.''', default=None
        )
    parser.add_argument(
        '--debug','-d',type=str,
        help='''(Optional) Provide an alternate location for the debug log and commands.''', default='logs/debuglog.txt'
        )
    args = parser.parse_args()
    if args.type in ('shortread','s','illumina'):
        args.type = 'shortread'
        fname = 'illumina_fastq/'
    elif args.type in ('longread','l','nanopore','ont'):
        args.type = 'longread'
        fname = 'ONT/'
    # append '/' to all args
    args.project,args.batch = [x+'/' if not x.endswith('/') else x for x in [args.project,args.batch]]
    if not os.path.isfile(args.input) or not os.path.isdir(args.project):
        print('Error locating provided files')
        quit(1)
    output_path = args.project + 'Sequence_data/' + fname + args.batch + 'raw_fastq/'
    # this is the path to the batch within the project folder
    # this should be something like /nfs/turbo/umms-esnitkin/Project_MDHHS_genomics/Sequence_data/illumina_fastq/2024-09-26_Plate1-to-Plate15/
    if not os.path.isdir('logs/'):
        subprocess.call(['mkdir','logs/'])
    with open(args.debug,'w') as debug:
        subprocess.call(['mkdir','-p',output_path],stdout=debug,stderr=debug)
        _ = debug.write(' '.join(['mkdir','-p',output_path]) + '\n')
        if args.qc is not None:
            if not args.qc.endswith('/'):
                args.qc+='/'
            if os.path.isdir(args.qc):
                # this moves the results folder from the QC pipeline and gives it a new name
                output_path_qc = args.project + 'Sequence_data/' + fname + args.batch + args.batch[:-1] + '_QC_results/'
                subprocess.call(['cp','-r',args.qc+'results/',output_path_qc])
                _ = debug.write(' '.join(['cp','-r',args.qc+'results/',output_path_qc]) + '\n')
            else:
                print('Warning: could not locate the provided directory of QC results.')
    move_and_rename(input_path=args.input,output_path=output_path,data_type=args.type,debuglog=args.debug)

if __name__ == "__main__":
    main()