read -p "Input uroborus detection result file:" sample 
if [ ! -e $sample ]; then
  echo "File not exists"
  exit
fi
awk '{print $1":"$2+1"|"$3"\t"$7}'     $sample > 1_circ.num && \
awk '$7>1{print $1":"$2+1"|"$3"\t"$7}' $sample > 2_circ.num 
