# return the path/sample name
def output(path,sample):            
	if len(path)<1:
		return sample
	else:
		if path[-1]=='/':
			return path+sample
		else:
			return path+'/'+sample


# obtain the tn barcode for each passed file
def exp_bc(file):
    bc_in=str(file).split('_')[-2]
    return bc_in
