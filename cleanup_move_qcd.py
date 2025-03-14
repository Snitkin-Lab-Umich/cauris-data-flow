import subprocess
import os
import argparse

# this script simply removes the large training files from the funannotate directory to allow easier transfer of files
# (expect this to take a while)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to the results directory generated after running a QCD pipeline. This will be in results/[batch_name]/
        Large intermediate files will be removed from this directory, and it will be moved to the output directory. 
        '''
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the new location where the QCD results should go. This should be an existing directory.
        '''
        )
    parser.add_argument(
        '--name','-n',type=str,
        help='''Provide the name for the new directory.
        '''
        )
    parser.add_argument(
        '--keep_training','-k',action='store_true',
        help='''If enabled, the training files will not be deleted.
        '''
        )
    args = parser.parse_args()
    if not os.path.isdir(args.input) or not os.path.isdir(args.output):
        print('Error locating provided directory')
        quit(1)
    args.input,args.output = [os.path.abspath(x)+'/' for x in [args.input,args.output]]
    if not args.keep_training:
        os.chdir(args.input)
        subprocess.call('rm -f funannotate/*/training/left.fq.gz',shell=True)
        subprocess.call('rm -f funannotate/*/training/right.fq.gz',shell=True)
        subprocess.call('rm -r -f funannotate/*/training/trimmomatic/',shell=True)
        subprocess.call('rm -r -f funannotate/*/training/trinity_gg/',shell=True)
    #print(' '.join(['mv',args.input,args.output + args.name]))
    subprocess.call(['mv',args.input,args.output + args.name])
    


if __name__ == "__main__":
    main()