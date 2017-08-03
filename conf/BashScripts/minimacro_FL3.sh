cd /Users/wilfeli/Documents
for wm in inf 1
do
	for learning in 0.05 0.1 0.5 0.9 0.95
	do
  		for seed in 2012 2013 2014 1 2 3 100 101 102 345 
  		do
     		input_file="minimacro.ini"
     		echo "N=1000" > $input_file
     		echo "seed=${seed}" >> $input_file
     		echo "F_QL_L=${learning}" >> $input_file
     		echo "H_QL_L=${learning}" >> $input_file
		echo "F_wm_LENGTH=${wm}" >> $input_file
		echo "H_wm_LENGTH=${wm}" >> $input_file
     		cat template_minimacro_FL.ini >> $input_file
     		./minimacro3 $input_file
  		done
	done
done
