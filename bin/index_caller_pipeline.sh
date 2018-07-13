script_folder=/Users/universe/Documents/2018_BG/index_caller/bin/
script_folder=/Users/gzx103/Documents/2018_spring/index_caller/bin/

##### get label info
cut -f1,2 atac_20cell.sig.18.txt | awk -F '_' -v OFS='\t' '{print $1"_"$2"_"$3"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$15"_"$20}' > position_index.txt
###### get signal matrix
cut -f3,5,6,7,8,9,10,11,15,20 atac_20cell.sig.18.txt > signal_mat.txt
###### paste label & signal_matrix
paste position_index.txt signal_mat.txt | awk -F '\t' -v OFS='\t' '{if ($2!="0_0_0_0_0_0_0_0_0_0") print $0}' > atac_20cell.sig.10.txt

	input_signal_mat = 'atac_20cell.sig.18.txt'
	output_signal_mat = 'data.adj.18.txt'
	index_count_lim = 100
	std_upper_lim = 3.0
	std_lower_lim = 0.1
	iteration_num_all = 100
	alpha_for_empty_index_set = 1.0
	random_seed = 2018




input_index_signal_mat=atac_20cell.sig.18.txt
ct_list=ct_list.18.txt

###### index-caller to recluster peaks
time python $script_folder'index_caller_mvn_bayes.py' -i $input_index_signal_mat -o $input_index_signal_mat'.ic.txt' -c 100 -u 3.0 -l 0.1 -t 100 -a 1.0 -r 2018

###### plot training curve (number of relabeled peaks for each iteration)
time Rscript $script_folder'plot_count_curve.R' $input_index_signal_mat'.ic.txt.change_num_all.txt' $input_index_signal_mat'.ic.txt.change_num_all.png' 2

###### relabel each index-set based on mean signal
time Rscript $script_folder'plot_pheatmap_index_caller.R' $input_index_signal_mat'.ic.txt.meansig.txt' $input_index_signal_mat'.ic.txt.meansig.png' $ct_list 2 red white F 0.0

###### replace old label by new label
time python $script_folder'oldlabel2newlabel.py' -i $input_index_signal_mat'.ic.txt' -l $input_index_signal_mat'.ic.txt.meansig.png.label2newlabel.txt' -c 2 -o $input_index_signal_mat'.newlabel.txt'

###### plot heatmap with new labels
sort -k2,2 $input_index_signal_mat'.newlabel.txt' > $input_index_signal_mat'.newlabel.sort.txt'
time Rscript $script_folder'plot_rect_sig.R' $input_index_signal_mat'.newlabel.sort.txt' $input_index_signal_mat'.newlabel.sort.png' $ct_list 3 dodgerblue white transparent F 0.01


sort -k2,2 $input_index_signal_mat'.ic.Nrelabel.txt' > $input_index_signal_mat'.ic.Nrelabel.sort.txt'
time Rscript $script_folder'plot_rect_sig.R' $input_index_signal_mat'.ic.Nrelabel.sort.txt' $input_index_signal_mat'.ic.Nrelabel.sort.png' $ct_list 3 dodgerblue white transparent F 0.01

time python $script_folder'relabel.py'
sort -k2,2 $input_index_signal_mat'.ic.Nrelabel1.txt' > $input_index_signal_mat'.ic.Nrelabel1.sort.txt'
time Rscript $script_folder'plot_rect_sig.R' $input_index_signal_mat'.ic.Nrelabel1.sort.txt' $input_index_signal_mat'.ic.Nrelabel1.sort.png' $ct_list 3 dodgerblue white transparent F 0.01
time Rscript $script_folder'plot_pheatmap_index_caller.R' $input_index_signal_mat'.ic.txt.meansig.txt' $input_index_signal_mat'.ic.txt.meansig.png' $ct_list 2 dodgerblue white F 0.0


sort -k2,2 $input_index_signal_mat'.ic.Nlabel.txt' > $input_index_signal_mat'.ic.Nlabel.sort.txt'
time Rscript $script_folder'plot_rect_sig.R' $input_index_signal_mat'.ic.Nlabel.sort.txt' $input_index_signal_mat'.ic.Nlabel.sort.png' $ct_list 3 dodgerblue white transparent F 0.01

###### plot heatmap of index-set mean signal
time Rscript $script_folder'plot_pheatmap_index_caller.R' $input_index_signal_mat'.ic.txt.meansig.txt' $input_index_signal_mat'.ic.txt.meansig.png' $ct_list 2 dodgerblue white F 0.0



###### get od signal matrix
time python $script_folder'get_od_index_set_meansig.py' -i $input_index_signal_mat -c 3 -o $input_index_signal_mat'.od.txt'

###### plot heatmap of index-set mean signal
time Rscript $script_folder'plot_pheatmap_index_caller.R' $input_index_signal_mat'.od.txt' $input_index_signal_mat'.od.png' $ct_list 2 dodgerblue white F 0.0

###### plot od heatmap
time Rscript $script_folder'plot_rect_sig.R' $input_index_signal_mat $input_index_signal_mat'.png' $ct_list 3 dodgerblue white transparent F 0.01


###### get peaks bed files
cell_type=$(tail -n+1 $ct_list | head -1 | cut -f1)
echo 1
echo $cell_type
cat $input_index_signal_mat'.newlabel.sort.txt' | awk -F '_' -v OFS='\t' -v colnum=$j '{print $1,$2,$3}' | awk -F '\t' -v OFS='\t' '{if ($4>0) print $1,$2,$3}' > $input_index_signal_mat'.'$cell_type'.new_pk.bed'
cat atac_20cell.sig.10.txt | awk -F '_' -v OFS='\t' -v colnum=$j '{print $1,$2,$3}' | awk -F '\t' -v OFS='\t' '{if ($4>0) print $1,$2,$3}' > $input_index_signal_mat'.'$cell_type'.old_pk.bed'
bedtools intersect -a $input_index_signal_mat'.'$cell_type'.new_pk.bed' -b $input_index_signal_mat'.'$cell_type'.old_pk.bed' -v > $input_index_signal_mat'.'$cell_type'.new_pk.spe.bed'
bedtools intersect -a $input_index_signal_mat'.'$cell_type'.old_pk.bed' -b $input_index_signal_mat'.'$cell_type'.new_pk.bed' -v > $input_index_signal_mat'.'$cell_type'.old_pk.spe.bed'

for i in {2..10}
do
	echo $i
	cell_type=$(tail -n+$i $ct_list | head -1 | cut -f1)
	j=$((i+2))
	echo $cell_type
	cat $input_index_signal_mat'.newlabel.sort.txt' | awk -F '_' -v OFS='\t' -v colnum=$j '{print $1,$2,$3,$colnum}' | awk -F '\t' -v OFS='\t' '{if ($5>0) print $1,$2,$3}' > $input_index_signal_mat'.'$cell_type'.new_pk.bed'
	cat $input_index_signal_mat | awk -F '_' -v OFS='\t' -v colnum=$j '{print $1,$2,$3,$colnum}' | awk -F '\t' -v OFS='\t' '{if ($5>0) print $1,$2,$3}' > $input_index_signal_mat'.'$cell_type'.old_pk.bed'
	bedtools intersect -a $input_index_signal_mat'.'$cell_type'.new_pk.bed' -b $input_index_signal_mat'.'$cell_type'.old_pk.bed' -v > $input_index_signal_mat'.'$cell_type'.new_pk.spe.bed'
	bedtools intersect -a $input_index_signal_mat'.'$cell_type'.old_pk.bed' -b $input_index_signal_mat'.'$cell_type'.new_pk.bed' -v > $input_index_signal_mat'.'$cell_type'.old_pk.spe.bed'
done


wc -l *.new_pk.bed > new_ct_pk_num.txt
wc -l *.old_pk.bed > old_ct_pk_num.txt


