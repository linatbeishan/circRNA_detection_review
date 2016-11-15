#threshold used: positive: 0-0, mixed-polyA: 0-0, background-simu(neg): 0.99-0.99, mixed-simu(pos+neg): 0.99-0.99, background-polyA: 0.97-0.97, Hela_plus:0.63-0.73, Hela_minus: 0.78-0.89, Hs68_minus: 0.78-0.75, Hs68_plus: 0.64-0.65, according to the cumalative distribution plot of posterior probability
perl glmReport-forCombine.pl Hela_RNaseR-plus_knife__circJuncProbs.txt 0.63 0.73 1>Hela_RNaseR-plus_1_circ.num 2> Hela_RNaseR-plus_2_circ.num
perl glmReport-forCombine.pl Hela_RNaseR-minus_knife__circJuncProbs.txt 0.78 0.89 1>Hela_RNaseR-minus_1_circ.num 2>Hela_RNaseR-minus_2_circ.num
perl glmReport-forCombine.pl Hs68_RNaseR-minus_knife__circJuncProbs.txt 0.78 0.75 1> Hs68_RNaseR-minus_1_circ.num 2> Hs68_RNaseR-minus_2_circ.num
perl glmReport-forCombine.pl Hs68_RNaseR-plus_knife__circJuncProbs.txt 0.64 0.65 1>Hs68_RNaseR-plus_1_circ.num 2>Hs68_RNaseR-plus_2_circ.num
perl glmReport-forCombine.pl positive_knife__circJuncProbs.txt 0 0 1>positive_1_circ.num 2>positive_2_circ.num
perl glmReport-forCombine.pl mixed-polyA_knife__circJuncProbs.txt 0 0 1>mixed-polyA_1_circ.num 2>mixed-polyA_2_circ.num
perl glmReport-forCombine.pl background-polyA_knife__circJuncProbs.txt 0.97 0.97 1>background-polyA_1_circ.num 2>background-polyA_2_circ.num
perl glmReport-forCombine.pl mixed-simu_knife__circJuncProbs.txt 0.99 0.99 mixed-simu_1_circ.num 2>mixed-simu_2_circ.num
perl glmReport-forCombine.pl background-simu_knife__circJuncProbs.txt 0.99 0.99 1>background-simu_1_circ.num 2>background-simu_2_circ.num
