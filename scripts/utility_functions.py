import glob

#adjust path with and without final slash
def output(path,sample):            
	if len(path)<1:
		return sample
	else:
		if path[-1]=='/':
			return path+sample
		else:
			return path+'/'+sample
#obtain sample name from sample folder (if a name is not specified)
def sample_name(string,path):   
    try:
        if len(string)<1:
            splt_path=path.split('/')
            if len(splt_path[-1])<1:
                sample_name=splt_path[-2]
                return sample_name
            else:
                sample_name=splt_path[-1]
                return sample_name
        else:
            sample_name=string
            return sample_name
    except IndexError:
        return 'scGET_sample'

# Generate the list of files with fastq.gz extension present in the input folder
def list_of_files(path,file):
	# input as a text file
	if len(file)>0:
		d={'1':[],'2':[],'3':[]}
		lof=[]
		fin=open(file).read().split('\n')[:-1]
		try:
			int(fin[0].strip().split(' ')[-1])
			for f in fin:
				nf=f.strip().split(' ')
				if len(path)<1:
					d[str(nf[-1])].append(path+nf[0])
				else:
					if path[-1]=='/':
						d[str(nf[-1])].append(path+nf[0])
					else:
						d[str(nf[-1])].append(path+'/'+nf[0])
		except ValueError:
			for f in fin:
				p=re.compile('R(1|2|3)')    # pattern indicating the read number (R1 or R2 or R3)
				nf=f.strip().split(' ')[0]      # file name
				r=p.search(f)[0][-1]        # r is the read number (i.e. p.search[0]=R1, while p.search[0][-1]=1)
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
    	# no input file (search in directory)
		if len(file)<1:
			if len(path)<1:
				return glob.glob('*fastq.gz')
			else:
				if path[-1]=='/':
					return glob.glob(path+'*fastq.gz')
				else:
					return glob.glob(path+'/'+'*fastq.gz')
		else:
    	# input: list with file names
			if str(file).startswith('['):
				list_input=[fin.strip() for fin in file]
				if len(path)<1:
					return list_input
				else:
					if path[-1]=='/':
						return [path+x for x in list_input]
					else:
						return [path+'/'+x for x in list_input]
    	# input: string with file names
			else:
				prelist=file.split(' ')
				list_input=[fin.strip() for fin in prelist]
				if len(path)<1:
					return list_input
				else:
					if path[-1]=='/':
						return [path+x for x in list_input]
					else:
						return [path+'/'+x for x in list_input]

# obtain the tn barcode for each passed file
def exp_bc(file):
    bc_in=str(file).split('_')[-2]
    return bc_in
