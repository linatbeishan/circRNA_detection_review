read -p "Input dcc detection result file:" sample 
if [ ! -e $sample ]; then
  echo "File not exists"
  exit
fi
awk 'NR>1{print $1":"$2"|"$3"\t"$4}' $sample >1_circ.num && \
awk 'NR>1&&$4>1{print $1":"$2"|"$3"\t"$4}' $sample >2_circ.num
