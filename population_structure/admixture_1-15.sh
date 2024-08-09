for iteration in {1..10}; do
  admixture --seed=$RANDOM --cv LK_intersect_norelated_wCro_nosexchrom.bed 15 > admx_LK_15\.$iteration.out
done

