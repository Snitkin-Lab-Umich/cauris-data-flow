import argparse
import os
import subprocess
import pandas as pd

def make_temp_reads_dir(input_dir, new_samples_file, temp_reads_dir):
    new_samples_df = pd.read_csv(new_samples_file)
    new_sample_cato_ids = new_samples_df['CaTO_ID'].tolist()
    # search through all source metadata files in input_dir
    # to avoid file collisions, make a dictionary of source file and destination file
    src_to_dest = {}
    for fname in os.listdir(input_dir):
        if fname.endswith('_sample_lookup.tsv'):
            source_metadata_file = os.path.join(input_dir, fname)
            source_metadata_df = pd.read_csv(source_metadata_file, sep='\t')
            for index, row in source_metadata_df.iterrows():
                fastq_path = row['fastq_path']
                sample_id = row['Sample_ID']
                cato_id = row['CaTO_ID']
                if cato_id not in new_sample_cato_ids:
                    continue
                r1_src = fastq_path.replace('_R*', '_R1')
                r2_src = fastq_path.replace('_R*', '_R2')
                r1_dest = os.path.join(temp_reads_dir, f'{sample_id}_R1.fastq.gz')
                r2_dest = os.path.join(temp_reads_dir, f'{sample_id}_R2.fastq.gz')
                if r1_src in src_to_dest or r2_src in src_to_dest or r1_dest in src_to_dest.values() or r2_dest in src_to_dest.values():
                    print(f'Error: file collision detected for {r1_src} or {r2_src}.')
                    quit(1)
                src_to_dest[r1_src] = r1_dest
                src_to_dest[r2_src] = r2_dest
    # copy files to temp_reads_dir
    for src in src_to_dest:
        dest = src_to_dest[src]
        if not os.path.isfile(src):
            print(f'Error: source file {src} does not exist.')
            quit(1)
        print(f'Copying {src} to {dest}')
        subprocess.run(['cp', src, dest])


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_dir','-i',type=str,
        help='''Provide a path to the directory of source metadata files.''',
        required=True
        )
    parser.add_argument(
        '--new_samples','-n',type=str,
        help='''Provide a path to the new sample file.''',
        required=True
        )
    parser.add_argument(
        '--temp_reads_dir','-t',type=str,
        help='''Provide a path to the directory where temporary reads will be stored.''',
        required=True
        )
    args = parser.parse_args()
    if not os.path.isdir(args.input_dir) or not os.path.isfile(args.new_samples):
        print(f'Could not locate directories or files')
        quit(1)
    if not os.path.isdir(args.temp_reads_dir):
        os.makedirs(args.temp_reads_dir)
    make_temp_reads_dir(args.input_dir, args.new_samples,args.temp_reads_dir)

if __name__ == "__main__":
    main()


