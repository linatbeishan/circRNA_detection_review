read -p "Input segemehl detection result file:" sample 
if [ ! -e $sample ]; then
  echo "File not exists"
  exit
fi
cat $sample |grep "C:"|grep -v "C:F" |grep -v "chrM" |awk -F ":" '{print $1,$2}' | sed 's/splits //g' |awk '{print $1":"$2"|"$3"\t"$4}'>1_circ.num && \
cat $sample |grep "C:"|grep -v "C:F" |grep -v "chrM" |awk -F ":" '{print $1,$2}' | sed 's/splits //g' |awk '$4>1{print $1":"$2"|"$3"\t"$4}' >2_circ.num
