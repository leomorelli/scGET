import pandas as pd
import os
import re
import gzip
import snakemake

#create_dirs() python function creates directories for each sample
def create_dirs(out_path,sample):
    sample_path = out_path+f'/{sample}'
    try:
        os.mkdir(out_path+f'/{sample}')
        os.mkdir(sample_path)
    except OSError:
        try:
            os.makedirs(sample_path)
            print("Successfully created the directory %s " % sample_path)
        except OSError:
            print ("Creation of the directory %s failed" % sample_path)
    else:
        print ("Successfully created the directory %s " % sample_path)
    return sample_path


#write_info() python function to create files for parallel multisample snakemake initialization

def write_info(out_path,info,sample):                                    
    with open(f'{out_path}/info_{sample}.txt','a') as myfile:
        for line in info:                                # info=information coming from samples, sample=sample name
            myfile.write(' '.join(map(str,line))+'\n')
    myfile.close()

#from a tex file to a list of files where files are classified according to their read number

def list_of_files(path,file):
        # input as a text file
        if len(file)>0:
                d={'1':[],'2':[],'3':[]}
                lof=[]
                fin=open(file).read().split('\n')[:-1]
                try:
                        int(fin[0].strip().split()[-1])
                        for f in fin:
                                nf=f.strip().split()
                                if len(path)<1:
                                        d[str(nf[-1])].append(path+nf[0])
                                else:   
                                        if path[-1]=='/':
                                                d[str(nf[-1])].append(path+nf[0])
                                        else:   
                                                d[str(nf[-1])].append(path+'/'+nf[0])
                except ValueError:
                        for f in fin:
                                nf=f.strip().split()[0]      # file name
                                r=re.search('(1|2|3)R_',nf[::-1])[0][0]        # (1|2|3)R_ pattern indicating the read number (R1 or R2 or R3), while r is the read number (i.e. p.search[0]=R1, while p.search[0][-1]=1)
                                if len(path)<1:
                                        d[str(r)].append(path+nf)
                                else:   
                                        if path[-1]=='/':
                                                d[str(r)].append(path+nf)
                                        else:   
                                                d[str(r)].append(path+'/'+nf)
                for i in d:
                        lof.append(d[i])
                return lof
        else:
            print('no input file provided')



def df_info(true_list):                                           # create a pd dataframe where the lane and the read of each file is clarified
        if str(true_list).startswith('[['):
                index=[]
                reads=[]
                for i in range(3):
                        for fin in true_list[i]:
                                index.append(fin)
                                reads.append(i+1)
                info=pd.DataFrame(reads,index=index,columns=['read'])
        else:
                info=pd.DataFrame(index=true_list,columns=['read'])     # index of dataframe represents file names. columns are lane and read
                for file in true_list:
                        fin=gzip.open(file,'rb').readline().decode().split()
                        info.loc[file,'read']=int(fin[-1][0])
        return info

def create_read_files(path,df):                                        # file names of file of the same reads are coupled in a text file
        for read in [1,2,3]:
                files=sorted(df[df.read==read].index)
                read_file=open(f'{path}/info_READ{read}.txt','w')              # text files are named 'merged_read?.txt'
                for line in files:
                        read_file.write(line)
                        read_file.write('\n')
                read_file.close()
#final function

def prep_input(input_file,out_path,input_path):
    df=pd.read_table(input_file,sep=' ', header=None,names=['file','read','name']) # the input file is converted into a pd.DataFrame where sample names, files names and read numbers are stored in specific columns
    samples=[x for x in df.loc[:,'name'].unique()]                               
    for smpl in samples:                                  # for each unique sample name a folder will be created containing information for sample analysis
        files_smpl=df.query('name == @smpl')
        text=[[x[0],x[1]] for x in files_smpl[['file', 'read']].values]
        smpl_path=create_dirs(out_path,smpl)        # create dirs
        write_info(smpl_path,text,smpl)             #write_info()
        lof=list_of_files(input_path,f'{smpl_path}/info_{smpl}.txt')
        df_smpl=df_info(lof)
        create_read_files(f'{smpl_path}',df_smpl)
