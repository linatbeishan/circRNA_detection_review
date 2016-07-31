read -p "Input ciri detection result file:" sample 
if [ ! -e $sample ]; then
  echo "File not exists"
  exit
fi
awk 'NR>1 && $5>0 {print $1"\t"$5}' $sample > 1_circ.num && \
awk 'NR>1 && $5>1 {print $1"\t"$5}' $sample > 2_circ.num
