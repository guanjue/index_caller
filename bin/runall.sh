input_dir='/Volumes/MAC_Data/data/labs/hardison_lab/vision/data_test/input_data/'
index_set_dir='/Volumes/MAC_Data/data/labs/hardison_lab/vision/data_test/index_set_matrix/'
index_set_sig_dir='/Volumes/MAC_Data/data/labs/hardison_lab/vision/data_test/index_set_sig_matrix/'

head $index_set_dir'celltype.index_set.sorted.txt'

head $input_dir'homerTable3.peaks.filtered.tpm.210k.txt'

head $index_set_dir'celltype.index.sorted.txt'

head $index_set_sig_dir'celltype.index_set.tpm.sorted.txt'

##################################
script_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/bin/'
index_caller_script_folder='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/index_caller/index_caller/bin/'

##################################
	###### initiate folders
	input_folder='input_data/'
	index_set_dir='index_set_matrix/'
	index_set_sig_dir='index_set_sig_matrix/'
	index_set_ideas_RE_dir='index_set_ideas_RE/'
	index_set_figure_dir='index_set_figure/'

	### mkdir index_set module folder
	if [ -d "$index_set_dir" ]; then  
    	rm -r $index_set_dir
	fi
	mkdir $index_set_dir

	### mkdir index_set module folder
	if [ -d "$index_set_sig_dir" ]; then  
    	rm -r $index_set_sig_dir
	fi
	mkdir $index_set_sig_dir

	### mkdir ideas state output folder
	if [ -d "$index_set_ideas_RE_dir" ]; then  
    	rm -r $index_set_ideas_RE_dir
	fi
	mkdir $index_set_ideas_RE_dir

	### mkdir figure folder
	if [ -d "$index_set_figure_dir" ]; then  
    	rm -r $index_set_figure_dir
	fi
	mkdir $index_set_figure_dir

##################################