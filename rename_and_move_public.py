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
            original_file = conv.iloc[i][2].strip()
            if not original_file.endswith(suffixes):
                print(f'Unexpected file format in {original_file}')
                quit(1)
            new_file = output_path + conv.iloc[i][1] + '.fastq.gz'
            with open(debuglog,'a') as debug:
                command = ['cp',original_file,new_file]
                subprocess.call(command,stdout=debug,stderr=debug)
                _ = debug.write(' '.join(command) + '\n')
        if data_type == 'shortread':
            search_str = conv.iloc[i][2].replace('*','_*').strip()
            # the conversion files have file names where a * immediately follows a number, which creates ambiguity between 1 and 10 and such
            # ex: ONT_Sequencing/Lojek_D5A_polished_results/Lojek_D5A_1*.fastq.gz and ONT_Sequencing/Lojek_D5A_polished_results/Lojek_D5A_10*.fastq.gz
            # this search_str line will need to be changed if the conversion files are formatted differently
            # note that .strip() is needed here, otherwise there can be issues with glob
            flist = sorted(glob.glob(search_str))
            if len(flist) != 2:
                print(f'Error while searching for files at path {conv.iloc[i][2]}')
                quit(1)
            original_file_r1,original_file_r2 = flist
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


def move_and_rename_public(input_path,output_path,sample_table,debuglog='logs/debuglog.txt'):
    # look in a QCD directory and move any of its fastq files from the input directory to the output directory
    samples_df = pd.read_csv(sample_table,sep='\t')
    if not os.path.isdir(output_path):
        subprocess.call(['mkdir','-p',output_path])
    with open(debuglog,'a') as debug:
        for sample in samples_df['sample_id']:
            search_path_r1 = os.path.join(input_path,f'{sample}_R1.fastq.gz')
            search_path_r2 = os.path.join(input_path,f'{sample}_R2.fastq.gz')
            if os.path.isfile(search_path_r1) and os.path.isfile(search_path_r2):
                cmd1 = ['cp',search_path_r1,os.path.join(output_path)]
                cmd2 = ['cp',search_path_r2,os.path.join(output_path)]
                _ = debug.write(' '.join(cmd1) + '\n')
                _ = debug.write(' '.join(cmd2) + '\n')
                subprocess.call(cmd1,stdout=debug,stderr=debug)
                subprocess.call(cmd2,stdout=debug,stderr=debug)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_path','-i',type=str,
        help='''Provide a path to the directory of raw reads used for this batch.'''
        )
    parser.add_argument(
        '--project','-p',type=str,
        help='''Provide a path to the project directory. This is the full file path that the QCD batch directories should go into.'''
        )
    parser.add_argument(
        '--batch','-b',type=str,
        help='''Provide a name for this batch of samples. A separate directory will be created with this name. Prefix this name with the date if possible.'''
        )
    parser.add_argument(
        '--samples','-s',type=str,
        help='''Provide a path to the samples.csv file used for this batch.'''
        )
    parser.add_argument(
        '--debug','-d',type=str,
        help='''(Optional) Provide an alternate location for the debug log and commands.''', default='logs/debuglog.txt'
        )
    args = parser.parse_args()
    if not os.path.isdir(args.input_path) or not os.path.isdir(args.project) or not os.path.isfile(args.samples):
        print('Error locating provided files')
        quit(1)
    output_path = os.path.join(args.project,args.batch,'raw_fastq/')
    # this is the path to the batch within the project folder
    # this should be something like /nfs/turbo/umms-esnitkin/Project_MDHHS_genomics/Sequence_data/illumina_fastq/2024-09-26_Plate1-to-Plate15/
    if not os.path.isdir('logs/'):
        subprocess.call(['mkdir','logs/'])
    with open(args.debug,'w') as debug:
        subprocess.call(['mkdir','-p',output_path],stdout=debug,stderr=debug)
        _ = debug.write(' '.join(['mkdir','-p',output_path]) + '\n')
    move_and_rename_public(input_path=args.input_path,output_path=output_path,sample_table=args.samples,debuglog=args.debug)

if __name__ == "__main__":
    main()