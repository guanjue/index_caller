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


data_old_label = read2d_array('data.adj.18.txt', str)

data_label2newlabel = read2d_array('qda_atac_wg.mvn_check.18.png.label2newlabel.txt', str)

data_label2newlabel_dcit = {}
for info in data_label2newlabel:
	data_label2newlabel_dcit[info[0]] = info[1]



data_new_label = []
for info in data_old_label:
	new_label = data_label2newlabel_dcit[info[1]]
	info[1] = new_label
	data_new_label.append(info)

data_new_label = np.array(data_new_label)

write2d_array(data_new_label, 'atac_20cell.sig.18.newlabel.txt')

