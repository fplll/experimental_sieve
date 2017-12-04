#!bin/bash
echo "Hello World!"

num_of_bases=5
dim_min=35
dim_max=60
ten=10

for ((dim=$dim_min; dim<=$dim_max; dim++)) 
do
	for i in {1..5} 
  do
		filename_out="Inputs/dim"$dim"_seed"$i
    #echo $filename_out
    #echo $((ten*dim))
    ./latticegen -randseed $i q $dim 1 $((ten*dim)) b >$filename_out
	done
done



