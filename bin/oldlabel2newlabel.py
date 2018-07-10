import os
import numpy as np
from scipy import stats

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return data0

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()


def oldlabel2newlabel(signal_mat_file, oldlabel2newlabel_mat_file, old_label_col, outputname):
	###### read signal with previous labels
	data_old_label = read2d_array(signal_mat_file, str)
	###### read old label to new label matrix
	data_oldlabel2newlabel = read2d_array(oldlabel2newlabel_mat_file, str)
	### get old label to new label dict
	data_oldlabel2newlabel_dcit = {}
	for info in data_oldlabel2newlabel:
		data_oldlabel2newlabel_dcit[info[0]] = info[1]

	###### replace old labels
	data_new_label = []
	for info in data_old_label:
		new_label = data_oldlabel2newlabel_dcit[info[old_label_col-1]]
		info[1] = new_label
		data_new_label.append(info)

	data_new_label = np.array(data_new_label)

	write2d_array(data_new_label, outputname)


############################################################################

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:l:c:o:")
	except getopt.GetoptError:
		print 'time python oldlabel2newlabel.py -i signal_mat_file -l oldlabel2newlabel_mat_file -c old_label_col_in_signal_mat -o outputname'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python oldlabel2newlabel.py -i signal_mat_file -l oldlabel2newlabel_mat_file -c old_label_col_in_signal_mat -o outputname'		
		elif opt=="-i":
			signal_mat_file=str(arg.strip())
		elif opt=="-l":
			oldlabel2newlabel_mat_file=str(arg.strip())
		elif opt=="-c":
			old_label_col=int(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())

	oldlabel2newlabel(signal_mat_file, oldlabel2newlabel_mat_file, old_label_col, outputname)

if __name__=="__main__":
	main(sys.argv[1:])

