script_folder=/Users/universe/Documents/2018_BG/index_caller/bin/
script_folder=/Users/gzx103/Documents/2018_spring/index_caller/bin/


input_index_signal_mat=atac_20cell.sig.18.txt
ct_list=ct_list.18.txt

cell_type=$(tail -n+1 $ct_list | head -1 | cut -f1)
echo 1
echo $cell_type
cat $input_index_signal_mat'.newlabel.sort.txt' | awk -F '_' -v OFS='\t' -v colnum=$j '{print $1,$2,$3}' | awk -F '\t' -v OFS='\t' '{if ($4>0) print $1,$2,$3}' > $input_index_signal_mat'.'$cell_type'.new_pk.bed'
cat atac_20cell.sig.10.txt | awk -F '_' -v OFS='\t' -v colnum=$j '{print $1,$2,$3}' | awk -F '\t' -v OFS='\t' '{if ($4>0) print $1,$2,$3}' > $input_index_signal_mat'.'$cell_type'.old_pk.bed'
bedtools intersect -a $input_index_signal_mat'.'$cell_type'.new_pk.bed' -b $input_index_signal_mat'.'$cell_type'.old_pk.bed' -v > $input_index_signal_mat'.'$cell_type'.new_pk.spe.bed'
bedtools intersect -a $input_index_signal_mat'.'$cell_type'.old_pk.bed' -b $input_index_signal_mat'.'$cell_type'.new_pk.bed' -v > $input_index_signal_mat'.'$cell_type'.old_pk.spe.bed'

for i in {2..18}
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


