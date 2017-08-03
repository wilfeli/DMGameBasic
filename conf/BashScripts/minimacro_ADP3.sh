cd /Users/wilfeli/Documents
for wm in 1 inf
do
	for grid in small big
	do
  		for seed in 2012 2013 2014 1 2 3 100 101 102 345 
  		do
     			input_file="minimacro.ini"
     			echo "N=1000" > $input_file
     			echo "seed=${seed}" >> $input_file
		       	echo "F_wm_LENGTH=${wm}" >> $input_file
			echo "H_wm_LENGTH=${wm}" >> $input_file
			echo "F_GRID=${grid}" >> $input_file
			echo "H_GRID=${grid}" >> $input_file
     			cat template_minimacro_ADP.ini >> $input_file
     			./minimacro3 $input_file
		done
	done
done
