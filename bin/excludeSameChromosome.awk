#!/usr/bin/awk -f

BEGIN {
  OFS="\t";
}
NR==FNR && $1 ~ /^>/ && $3 ~/^chromosome/ {
  #gsub("^>","",$1)
  split($3,location,":");
  idmap[$1]=location[3];
}
#PROCESSING PAIRS, print all across datasets, but skip if same dataset and same block (chromosome)
# NR!=FNR{print $1,$2}
NR!=FNR && $1 != $2 && (tag1 != tag2 || ($1 in idmap && $2 in idmap && idmap[$1] != idmap[$2])) {
  print;
} 