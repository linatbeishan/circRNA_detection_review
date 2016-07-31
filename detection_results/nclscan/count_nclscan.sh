read -p "Input nclscan detection result file:" sample 
if [ ! -e $sample ]; then
  echo "File not exists"
  exit
fi
#intragenic
awk '$9==1 {if($3=="+") print $1":"$5"|"$2"\t"$11;  else print $1":"$2"|"$5"\t"$11}' $sample > 1_circ.num && \
#intragenic & numJuncs > 1
awk '$9==1 && $11>1{if($3=="+") print $1":"$5"|"$2"\t"$11; else print $1":"$2"|"$5"\t"$11}' $sample > 2_circ.num
