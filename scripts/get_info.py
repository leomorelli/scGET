import gzip
import yaml
import pandas as pd
import importlib.util

spec = importlib.util.spec_from_file_location("utility_functions", "scripts/utility_functions.py")
utilities = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utilities)

yaml_file=open('config.yaml')
config=yaml.load(yaml_file,Loader=yaml.FullLoader)
output_path=utilities.output(config['output_path'])


def filter_files(files_list):                                     # filter only READ1,2 and 3 files from the folder with input files 
	true_list=list()
	for file in files_list:
		fin=gzip.open(file,'rb').readline().decode().split(' ')   # list of words of the first line of the file
		if len(fin)!=2:                                           # sequencing files contain two words in the first line
			continue
		else:
			try:
				if int(fin[-1][0]) in [1,2,3]:                    # the last word must begin with a number between 1,2,3 (read number)
					true_list.append(file)
			except ValueError:
				continue
	return true_list

def df_info(true_list):                                           # create a pd dataframe where the lane and the read of each file is clarified
	info=pd.DataFrame(index=true_list,columns=['read'])    # index of dataframe represents file names. columns are lane and read
	for file in true_list:
		fin=gzip.open(file,'rb').readline().decode().split(' ')
		info.loc[file,'read']=int(fin[-1][0])
	return info

def create_read_files(df,path):                                        # file names of file of the same reads are coupled in a text file
	for read in [1,2,3]:
		files=info[info.read==read].index
		read_file=open(f'{path}/info_READ{read}.txt','w')              # text files are named 'merged_read?.txt'
		for line in files:
			read_file.write(line)
			read_file.write('\n')
		read_file.close()

info=df_info(filter_files(snakemake.input))


#info.to_csv('info_read.csv')

create_read_files(info,output_path)