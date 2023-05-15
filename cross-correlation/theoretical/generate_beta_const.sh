beta=1.0
#c=0

for z0 in 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.22;
do
	for lbda in 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8;
	do
		#(( c ++ ))
		./bin/pk2ctg -pk Cl_Pk/pk_3dmatter.dat -wl 0.6847 -wm 0.3153 -h 0.6736 -lmax 128 -z0 $z0 -beta $beta -lbda $lbda -bg 1 -out funcs_beta_$beta/z0_$z0.lbda_$lbda.dat
		echo "z0=$z0, beta=$beta e lbda=$lbda"
		#echo "$c combinações concluidas (total=100)"
	done
done
