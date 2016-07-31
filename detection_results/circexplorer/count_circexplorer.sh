read -p "Input circexplorer detection result file:" sample 
if [ ! -e $sample ]; then
  echo "File not exists"
  exit
fi
awk '$13>0{print $1":"$2+1"|"$3"\t"$13}' $sample >1_circ.num && \
awk '$13>1{print $1":"$2+1"|"$3"\t"$13}' $sample >2_circ.num
