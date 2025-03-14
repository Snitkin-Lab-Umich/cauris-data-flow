import subprocess
import os
import argparse

# this script moves the hybrid assembly directory to a new location
# it also copies the main outputs to the central directories
# currently, these are: [polypolish,funannotate_update,multiqc]

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
    os.chdir(args.input)
    if not args.keep_training:
        remove_training_files()
    #print(' '.join(['mv',args.input,args.output + args.name]))
    #subprocess.call(['mv',args.input,args.output + args.name])
    move_and_copy(results_dir=args.input,output_dir=args.output,output_name=args.name)
    
def remove_training_files():
    subprocess.call('rm -f funannotate/*/training/left.fq.gz',shell=True)
    subprocess.call('rm -f funannotate/*/training/right.fq.gz',shell=True)
    subprocess.call('rm -r -f funannotate/*/training/trimmomatic/',shell=True)
    subprocess.call('rm -r -f funannotate/*/training/trinity_gg/',shell=True)

def move_and_copy(results_dir,output_dir,output_name):
    new_results_dir = output_dir + output_name + '/'
    # move entire directory to new location
    command1 = ['mv',results_dir,new_results_dir]
    subprocess.call(command1)
    # copy the contents of each folder to the main directories
    # prompt before adding duplicate data
    command2 = f'cp -rib {new_results_dir}polypolish/*/ {output_dir}polypolish/'
    subprocess.run(command2,shell=True)
    # this removes the large and unnecessary files present in the polypolish output
    command22 = f'rm {output_dir}polypolish/*/*.sam'
    subprocess.run(command22,shell=True)
    command3 = f'cp -rib {new_results_dir}quast/*/ {output_dir}quast/'
    subprocess.run(command3,shell=True)
    command4 = f'cp -rib {new_results_dir}funannotate/*/ {output_dir}funannotate/'
    subprocess.run(command4,shell=True)
    print(command1)
    print(command2)
    print(command22)
    print(command3)
    print(command4)

if __name__ == "__main__":
    main()