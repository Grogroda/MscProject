z0=0.054
c=0

for beta in 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0;
do
	for lbda in 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0;
	do
		(( c ++ ))
		./bin/pk2ctg -pk Cl_Pk/pk_3dmatter.dat -wl 0.6847 -wm 0.3153 -h 0.6736 -lmax 128 -z0 $z0 -beta $beta -lbda $lbda -bg 1 -out funcs_z0_0.054/beta_$beta.lbda_$lbda.dat
		echo "$c combinações concluidas (total=100)"
	done
done
