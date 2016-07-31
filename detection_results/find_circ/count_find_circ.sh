read -p "Input find_circ detection result file:" sample 
if [ ! -e $sample ]; then
  echo "File not exists"
  exit
fi
grep -v "chrM" $sample |awk '$5>1{print $1":"$2+1"|"$3"\t"$5}' > 2_circ.num && \
grep -v "chrM" $sample |awk '$5>0{print $1":"$2+1"|"$3"\t"$5}' > 1_circ.num
