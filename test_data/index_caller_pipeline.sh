script_folder=/Users/universe/Documents/2018_BG/index_caller/bin/

##### get label info
cut -f1,2 atac_20cell.sig.18.txt | awk -F '_' -v OFS='\t' '{print $1"_"$2"_"$3"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$15"_"$20}' > position_index.txt
###### get signal matrix
cut -f3,5,6,7,8,9,10,11,15,20 atac_20cell.sig.18.txt > signal_mat.txt
###### paste label & signal_matrix
paste position_index.txt signal_mat.txt | awk '{if ($2!="0_0_0_0_0_0_0_0_0_0") print $0}' > atac_20cell.sig.10.txt

	input_signal_mat = 'atac_20cell.sig.10.txt'
	output_signal_mat = 'data.adj.10.txt'
	index_count_lim = 100
	std_upper_lim = 3.0
	std_lower_lim = 0.1
	iteration_num_all = 100
	alpha_for_empty_index_set = 1.0
	random_seed = 2018

###### index-caller to recluster peaks
time python $script_folder'index_caller_mvn_bayes.py' -i atac_20cell.sig.10.txt -o data.adj.10.txt -c 100 -u 3.0 -l 0.1 -t 100 -a 1.0 -r 2018

###### plot training curve (number of relabeled peaks for each iteration)
time Rscript $script_folder'plot_count_curve.R' data.adj.10.txt.change_num_all.txt data.adj.10.txt.change_num_all.png 2

###### relabel each index-set based on mean signal
time Rscript $script_folder'plot_pheatmap_index_caller.R' data.adj.10.txt.meansig.txt data.adj.10.txt.meansig.png ct_list.10.txt 2 red white F 0.0

###### replace old label by new label
time python $script_folder'oldlabel2newlabel.py' -i data.adj.10.txt -l data.adj.10.txt.meansig.png.label2newlabel.txt -c 2 -o atac_20cell.sig.10.newlabel.txt

###### plot heatmap with new labels
sort -k2,2 atac_20cell.sig.10.newlabel.txt > atac_20cell.sig.10.newlabel.reorder.txt
time Rscript $script_folder'plot_rect_sig.R' atac_20cell.sig.10.newlabel.reorder.txt atac_20cell.sig.10.newlabel.reorder.png ct_list.10.colnames.txt 3 dodgerblue white transparent F 0.01

###### plot heatmap of index-set mean signal
time Rscript $script_folder'plot_pheatmap_index_caller.R' data.adj.10.txt.meansig.txt data.adj.10.txt.meansig.png ct_list.10.colnames.txt 2 dodgerblue white F 0.0



###### plot heatmap of index-set mean signal
time Rscript $script_folder'plot_pheatmap_index_caller.R' qda_atac_wg.mvn_check.10.od.txt qda_atac_wg.mvn_check.10.od.png ct_list.10.txt 2 dodgerblue white F 0.0

###### get od signal matrix
time python $script_folder'get_od_index_set_meansig.py' -i atac_20cell.sig.10.txt -c 3 -o qda_atac_wg.mvn_check.10.od.txt

###### plot od heatmap
time Rscript $script_folder'plot_rect_sig.R' atac_20cell.sig.10.txt atac_20cell.sig.10.png ct_list.10.colnames.txt 3 dodgerblue white transparent F 0.01




