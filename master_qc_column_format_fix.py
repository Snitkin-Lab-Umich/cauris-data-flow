import pandas as pd
import subprocess
import os
import argparse

def hamming_distance(str1, str2):
    # calculate the hamming distance between two strings
    if len(str1) != len(str2):
        print("Strings are not of equal length")
        return(len(str1) + len(str2))
    return sum(el1 != el2 for el1, el2 in zip(str1, str2))

def column_format_fix(input, output):
    # make a backup of the output file
    if os.path.isfile(output):
        backup_file = output.replace('.csv','_backup.csv')
        subprocess.run(['cp',output,backup_file])
        print(f'Backup of {output} created at {backup_file}')
    df_in = pd.read_csv(input, sep='\t')
    df_out = pd.read_csv(output, sep='\t')
    # remove the collection column from the output dataframe if it exists
    if 'collection' in df_out.columns:
        df_out = df_out.drop(columns=['collection'])
    # move the data_type column the position immediately after auriclass_clade
    if 'data_type' in df_out.columns:
        col_list = df_out.columns.tolist()
        col_list.remove('data_type')
        auriclass_index = col_list.index('auriclass_clade')
        col_list.insert(auriclass_index + 1, 'data_type')
        df_out = df_out.reindex(columns=col_list)
    # compare each column to make sure they have a similar hamming distance
    if len(df_in.columns) != len(df_out.columns):
        print(f'Error: The number of columns in {input} ({len(df_in.columns)}) does not match the number of columns in {output} ({len(df_out.columns)}).')
        quit(1)
    for col_in, col_out in zip(df_in.columns, df_out.columns):
        if hamming_distance(col_in, col_out) > 4:
            if 'contig' in col_in and 'contig' in col_out:
                continue
            if 'GC' in col_in and 'GC' in col_out:
                continue
            print(f'Error: Differing column names found: {input}: {col_in}, {output}: {col_out}')
            quit(1)
    # if everything looks close, simply replace the output's columns with the input's columns
    df_out.columns = df_in.columns
    df_out.to_csv(output, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to an input qc file with the correct column names.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a path an output qc file to modify.''',
        required=True
        )
    args = parser.parse_args()
    column_format_fix(input=args.input, output=args.output)

if __name__ == "__main__":
    main()







