lbda=3.4
c=0

for z0 in 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.22;
do
	for beta in 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4;
	do
		(( c ++ ))
		./bin/pk2ctg -pk Cl_Pk/pk_3dmatter.dat -wl 0.6847 -wm 0.3153 -h 0.6736 -lmax 128 -z0 $z0 -beta $beta -lbda $lbda -bg 1 -out funcs_lbda_$lbda/z0_$z0.beta_$beta.dat
		echo "$c combinações concluidas (total=63)"
	done
done

