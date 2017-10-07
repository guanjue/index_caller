time python reads_count_tpm.py -i reads_count_matrix_5end_whole_genome.add_id.txt -o reads_count_matrix_5end_whole_genome_tpm.txt

time python cell_type_sort.py -i reads_count_matrix_5end_whole_genome_tpm.txt -b bam_file.txt -r homerTable3.peaks.filtered.txt -o homerTable3.peaks.filtered.tpm.txt

time python merge_cell_type_data.py -i homerTable3.peaks.filtered.tpm.txt -m sample2celltype.txt -n 2 -o celltype.tpm.txt



cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$2}' > lsk.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$3}' > cmp.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$4}' > gmp.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$5}' > mep.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$6}' > cfue.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$7}' > ery.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$8}' > cfumk.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$9}' > meg.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$10}' > mono.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$11}' > neu.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$12}' > b.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$13}' > nk.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$14}' > tcd4.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$15}' > tcd8.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$16}' > g1e.wg.tpm.txt
cat celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$17}' > er4.wg.tpm.txt







cat DNA_regin_210k_indexsort_onlyinterval.txt | awk -F '\t' -v OFS='\t' '{print $2,$3,$4}' | sort -k1,1 -k2,2n > 210k_sorted.bed

paste wg.bed celltype.tpm.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3, $5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$15"_"$16"_"$17"_"$18"_"$19"_"$20}' > celltype.tpm.bed

bedtools intersect -a celltype.tpm.bed -b 210k_sorted.bed -v > celltype.tpm.NOcRE.bed

cat celltype.tpm.NOcRE.bed | awk -F '_' -v OFS='\t' '{print $1,$2, $4,$3, $5,$6, $15,$16, $7,$8,$9,$10,$11,$12,$13,$14}' > celltype.tpm.NOcRE.txt
tail -n+2 celltype.tpm.bed | awk -F '_' -v OFS='\t' '{print $1,$2, $4,$3, $5,$6, $15,$16, $7,$8,$9,$10,$11,$12,$13,$14}' > celltype.tpm.sorted.txt



