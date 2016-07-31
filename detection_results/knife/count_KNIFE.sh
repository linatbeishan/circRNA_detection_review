#!/bin/bash

read -p "Input KNIFE detection result file:" sample 

case=1
if [ ! -e $sample ]; then
  echo "File not exists"
  exit
elif [[ $sample =~ .*RNaseR-minus_knife__circJuncProbs.* ]]; then
  program="glmReport.pl"
  th=0.78
elif [[ $sample =~ .*RNaseR-plus_knife__circJuncProbs.* ]]; then
  program="glmReport.pl"
  th=0.62
elif [[ $sample =~ ^polyA_knife__circJuncProbs.* ]]; then
  program="glmReport.pl"
  th=0.90
elif [[ $sample =~ ^mixed_knife_report.* ]] || [[ $sample =~ ^positive_knife_report.* ]]; then
  program="report.pl";
  case=0
else
  echo "Input file error, please check the usage of: glmReport.pl and report.pl"
  exit
fi

if [ $case == 1 ]; then
  echo -e "perl $program $sample $th 1>1_circ.num 2>2_circ.num";
  perl $program $sample $th 1>1_circ.num 2>2_circ.num;
else
  echo -e "perl $program $sample 1>1_circ.num 2>2_circ.num";
  perl $program $sample 1>1_circ.num 2>2_circ.num;
fi
