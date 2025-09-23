import pandas as pd
import subprocess
import os
import argparse

def move_between_dirs(source_dir,dest_dir,source_qc_file,dest_qc_file,temp_qc_file,debuglog = 'logs/debug.txt'):
    source_dir,dest_dir = [x+'/' if not x.endswith('/') else x for x in [source_dir,dest_dir]]
    # for each directory in source_dir/funannotate, check if it exists in dest_dir/funannotate
    # if it does not AND if it is in the qc_file as PASS, move all of its directories over and add it to the temp_qc_file
    source_master_qc = pd.read_csv(source_qc_file, sep = '\t')
    dest_master_qc = pd.read_csv(dest_qc_file, sep = '\t')
    isolates_to_save = []
    with open(debuglog, 'a') as debug:
        for isolate_name in os.listdir(source_dir + 'funannotate/'):
            if os.path.isdir(source_dir + 'funannotate/' + isolate_name):
                source_qc = check_pass(source_master_qc,isolate_name)
                dest_qc = check_pass(dest_master_qc,isolate_name)
                #print([source_qc,dest_qc])
                if source_qc != 'present_pass':
                    print(f'Error with {isolate_name} - present in assembly directory but not marked as PASS in the corresponding source QC file')
                    quit(1)
                if dest_qc == 'absent' or dest_qc == 'present_fail':
                    if not os.path.isdir(dest_dir + 'funannotate/' + isolate_name):
                        cmdlist = []
                        for subdir in ['funannotate/','busco/','quast/','spades/']:
                            cmdlist.append(['cp','-r',source_dir + subdir + isolate_name,dest_dir + subdir + isolate_name])
                        for cmd in cmdlist:
                            subprocess.call(cmd,stdout=debug,stderr=debug)
                            _ = debug.write(' '.join(cmd) + '\n')
                        # add this row to the list of rows to save
                        isolates_to_save.append(isolate_name)
                    else:
                        print(f'Error with {isolate_name} - present in assembly directory but not marked as PASS in the corresponding destination QC file')
                        quit(1)
    # now, subset the source_master_qc to just the isolates_to_save and save as a temp_qc_file
    if len(isolates_to_save) > 0:
        temp_df = source_master_qc[source_master_qc['Sample'].isin(isolates_to_save)]
        temp_df.to_csv(temp_qc_file, sep = '\t', index = False)


def check_pass(qc_df,sample):
    if sample in qc_df['Sample'].values:
        row = qc_df[qc_df['Sample'] == sample]
        if row['QC_EVALUATION'].values[0] == 'PASS':
            return('present_pass')
        else:
            return('present_fail')
    else:
        return('absent')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--source_dir','-sd',type=str,
        help='''Provide a path to the source directory to copy from. This should be an assembly/ directory containing 
        funannotate/, busco/, quast/, and spades/ subdirectories.
        ''')
    parser.add_argument(
        '--dest_dir','-dd',type=str,
        help='''Provide a path to the source directory to copy from. This should be an assembly/ directory containing 
        funannotate/, busco/, quast/, and spades/ subdirectories.
        ''')
    parser.add_argument(
        '--source_qc','-sqc',type=str,
        help='''Provide a path to the master QC summary file for the source directory.
        ''')
    parser.add_argument(
        '--dest_qc','-dqc',type=str,
        help='''Provide a path to the master QC summary file for the destination directory.
        ''')
    parser.add_argument(
        '--temp_qc','-tqc',type=str,
        help='''Provide a path to a temporary QC summary file that contains only the isolates that were moved.
        ''')
    args = parser.parse_args()
    if not os.path.isdir('logs/'):
        os.mkdir('logs/')
    move_between_dirs(
        source_dir = args.source_dir, dest_dir = args.dest_dir, source_qc_file = args.source_qc, dest_qc_file = args.dest_qc, 
        temp_qc_file = args.temp_qc,debuglog = 'logs/debug.txt')

if __name__ == "__main__":
    main()







