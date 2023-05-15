echo "Directory of the functions:"

read dir ;

for file in $dir/* ;
do
	mv $file $dir/temporary.dat
	#ls $dir/
	sort -n -k 1 $dir/temporary.dat > $file
	#ls $dir/
	rm $dir/temporary.dat
	#ls $dir/
done;
