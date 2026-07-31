import argparse
import os
import pandas as pd

def make_fastq_dict(new_df, reads_dir, r1_suffix, r2_suffix):
    sid_to_fastq = {}
    names_used = set()
    for index, row in new_df.iterrows():
        sid = row['Sequencing_ID']
        fastq_names = []
        for fname in os.listdir(reads_dir):
            if fname.startswith(sid) and (fname.endswith(r1_suffix) or fname.endswith(r2_suffix)):
                fastq_names.append(fname)
        if len(fastq_names) == 2 and sum(fname.endswith(r1_suffix) for fname in fastq_names) == 1 and sum(fname.endswith(r2_suffix) for fname in fastq_names) == 1:
            if any(fname in names_used for fname in fastq_names):
                print(f'Error: fastq files for {sid} have already been used. File list: {fastq_names}')
                quit(1)
            names_used.update(fastq_names)
            sid_to_fastq[sid] = fastq_names
        else:
            if len(fastq_names) == 0:
                fastq_names = 'No files found'
            print(f'Error locating fastq files for {sid}. File list: {fastq_names}')
            quit(1)
    return(sid_to_fastq)


def process_samples(input_csv, reads_dir, output_dir, seqtype, alt_suffix):
    new_df = pd.read_csv(input_csv)
    if new_df.columns.tolist() != ['CaTO_ID', 'Sequencing_ID', 'Source']:
        print(f'The columns of the input csv must be CaTO_ID, Sequencing_ID, and Source')
        quit(1)
    if alt_suffix is None:
        r1_suffix = '_R1.fastq.gz'
        r2_suffix = '_R2.fastq.gz'
        gen_suffix = '_R*.fastq.gz'
    else:
        r1_suffix = alt_suffix
        r2_suffix = alt_suffix.replace('R1','R2')
        gen_suffix = alt_suffix.replace('R1','R*')
    # to avoid file collisions, make a dictionary for each Sequencing_ID and its corresponding R1 and R2 fastq.gz files
    sid_to_fastq = make_fastq_dict(new_df, reads_dir, r1_suffix, r2_suffix)
    # group by source and assign each sample to its corresponding source metadata file
    for source, new_df_source in new_df.groupby('Source'):
        # make sure the source metadata file is ok
        source_metadata_file = os.path.join(output_dir, f'{source}_{seqtype}_sample_lookup.tsv')
        if not os.path.isfile(source_metadata_file):
            print(f'Error: could not locate source metadata file {source_metadata_file}')
            quit(1)
        source_metadata_df = pd.read_csv(source_metadata_file,sep='\t')
        if source_metadata_df.columns.tolist() != ['CaTO_ID', 'Sample_ID', 'fastq_path']:
            print(f'The columns of the source metadata csv must be CaTO_ID, Sample_ID, and fastq_path')
            print(f'Columns found: {source_metadata_df.columns.tolist()}')
            quit(1)
        source_metadata_df = pd.read_csv(source_metadata_file, sep = '\t')
        # Sample_ID should be in the format [source]_Caur_[number]
        # determine highest number present in Sample_ID column
        sample_numbers = [int(x.split('_Caur_')[-1]) for x in source_metadata_df['Sample_ID'].tolist()]
        starting_index = max(sample_numbers)
        # append each new sample to the source metadata file with a unique Sample_ID and its corresponding fastq_path
        current_index = starting_index + 1
        for index, row in new_df_source.iterrows():
            cato_id = row['CaTO_ID']
            sid = row['Sequencing_ID']
            fastq_path = sid_to_fastq[sid][0]
            fastq_path = os.path.join(reads_dir, fastq_path)
            fastq_path = fastq_path.replace(r1_suffix, gen_suffix)
            fastq_path = fastq_path.replace(r2_suffix, gen_suffix)
            while current_index in sample_numbers:
                current_index += 1
            sample_id = f'{source}_Caur_{current_index}'
            sample_numbers.append(current_index)
            new_row = pd.DataFrame([[cato_id, sample_id, fastq_path]], columns=['CaTO_ID', 'Sample_ID', 'fastq_path'])
            source_metadata_df = pd.concat([source_metadata_df, new_row], ignore_index=True)
            current_index += 1
        # write the updated source metadata file
        source_metadata_df.to_csv(source_metadata_file, sep = '\t', index=False)
        print(f'Updated source metadata file {source_metadata_file}')
        print(f'Added samples: {new_df_source["Sequencing_ID"].tolist()}')


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to a csv file with the samples to add. The csv should have columns for CaTO_ID, Sequencing_ID, and Source.''',
        required=True
        )
    parser.add_argument(
        '--output_dir','-o',type=str,
        help='''Provide a path to the output directory. This should contain metadata files for each source.''',
        required=True
        )
    parser.add_argument(
        '--reads_dir','-r',type=str,
        help='''Provide a path to the directory containing the raw reads. Files to process must end in _R1.fastq.gz or _R2.fastq.gz. Sample names must start with the Sequencing_ID column in the input csv.''',
        required=True
        )
    parser.add_argument(
        '--seqtype','-s',type=str,choices=['shortread','longread'],
        help='''Provide the sequencing type. Must be either 'shortread' or 'longread'.''',
        required=True
        )
    parser.add_argument(
        '--alt_suffix','-a',type=str,
        help='''If the fastq files do not end in _R1.fastq.gz or _R2.fastq.gz, provide the alternative suffix here. Provide this suffix as it is present in the R1 reads.
        Example: _R1_001.fastq.gz''',
        required=False,default=None
        )
    args = parser.parse_args()
    if not os.path.isdir(args.reads_dir) or not os.path.isfile(args.input) or not os.path.isdir(args.output_dir):
        print(f'Could not locate directories or files')
        quit(1)
    process_samples(args.input, args.reads_dir, args.output_dir, args.seqtype, args.alt_suffix)

if __name__ == "__main__":
    main()


