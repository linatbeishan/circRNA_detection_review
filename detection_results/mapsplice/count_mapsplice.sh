read -p "Input mapsplice detection result file:" sample 
if [ ! -e $sample ]; then
  echo "File not exists"
  exit
fi
cut -d"~" -f2 $sample |awk '$5>0{ if($6=="++") print $1":"$3"|"$2"\t"$5; else print $1":"$2"|"$3"\t"$5; }' > 1_circ.num && \
cut -d"~" -f2 $sample |awk '$5>1{ if($6=="++") print $1":"$3"|"$2"\t"$5; else print $1":"$2"|"$3"\t"$5; }' > 2_circ.num
