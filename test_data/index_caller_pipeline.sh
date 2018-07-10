cut -f1,2 atac_20cell.sig.18.txt | awk -F '_' -v OFS='\t' '{print $1"_"$2"_"$3"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$15"_"$20}' > position_index.txt

cut -f3,5,6,7,8,9,10,11,15,20 atac_20cell.sig.18.txt > signal_mat.txt


paste position_index.txt signal_mat.txt | awk '{if ($2!="0_0_0_0_0_0_0_0_0_0") print $0}' > atac_20cell.sig.10.txt


time python /Users/gzx103/Documents/2018_spring/index_caller/bin/index_caller_mvn_bg.py


time Rscript /Users/gzx103/Documents/2018_spring/index_caller/bin/plot_count_curve.R change_num_all.10.txt change_num_all.10.counts.png 2


time Rscript /Users/gzx103/Documents/2018_spring/index_caller/bin/plot_pheatmap_index_caller.R qda_atac_wg.mvn_check.10.txt qda_atac_wg.mvn_check.10.png ct_list.10.txt 2 red white F 0.0


time python /Users/gzx103/Documents/2018_spring/index_caller/bin/index_caller_mvn_od_mean_mat.py


time Rscript /Users/gzx103/Documents/2018_spring/index_caller/bin/plot_pheatmap_index_caller.R qda_atac_wg.mvn_check.10.od.txt qda_atac_wg.mvn_check.10.od.png ct_list.10.txt 2 red white F 0.0


time Rscript plot_rect_sig.R atac_20cell.sig.od.txt atac_18cell.sig.od.png ct_list.18.colnames.txt 3 red white transparent F 0.01

time Rscript plot_rect_sig.R atac_20cell.sig.18.newlabel.reorder.txt atac_20cell.sig.18.newlabel.reorder.png ct_list.18.colnames.txt 3 dodgerblue white transparent F 0.01
