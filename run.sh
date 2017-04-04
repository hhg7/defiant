cd /home/con/DNA_Methylation/Public_data/radmeth/rrbs/1/1
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/rrbs/1/2
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/rrbs/1/3
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/rrbs/1/4
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/rrbs/2/1
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/rrbs/2/2
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/rrbs/2/3
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/rrbs/2/4
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/wgbs/1/1
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/wgbs/1/2
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/wgbs/1/3
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/wgbs/1/4
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/wgbs/2/1
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/wgbs/2/2
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/wgbs/2/3
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
cd /home/con/DNA_Methylation/Public_data/radmeth/wgbs/2/4
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/wand -factor case /home/con/DNA_Methylation/Public_data/radmeth/design_matrix.tsv proportion_table.tsv > cpgs.bed 2> wand.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/adjust -bins 1:100:1 cpgs.bed > cpgs.adjusted.bed 2> adjust.stderr
/usr/bin/time awk '$5 < 0.01 "{ print  ; $}"' cpgs.adjusted.bed > dm_cpgs.bed 2> awk.stderr
/usr/bin/time /home/con/DNA_Methylation/Public_data/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed 2> dmrs.stderr
