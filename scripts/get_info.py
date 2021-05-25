import gzip
import pandas as pd



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
