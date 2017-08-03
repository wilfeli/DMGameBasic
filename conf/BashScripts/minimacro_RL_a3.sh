cd /Users/wilfeli/Documents
for fopt in mRE
do
	for hopt in QL EO_CS EO_ADP
	do
		for seed in 2012 2013 2014 1 2 3 100 101 102 345 
		do
			input_file="minimacro.ini"
			echo "N=1000" > $input_file
			echo "seed=${seed}" >> $input_file
			echo "F_opt_TYPE=${fopt}" >> $input_file
			echo "H_opt_TYPE=${hopt}" >> $input_file	
			cat template_minimacro_RL_a.ini >> $input_file
			./minimacro3 $input_file
		done
	done
done
