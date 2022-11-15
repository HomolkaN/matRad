#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH --job-name=topas_testrun$1      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=28            # Number of CPU cores per task
#SBATCH --output=topas_testrun_%j.log     # Standard output and error log
pwd; hostname; date

export TOPAS_G4_DATA_DIR=/opt/topas/3.8.1/G4Data


counter=0
for f in $(find $PWD/$1 -maxdepth 3 -type d | sort -n)
do
	res=($f/*run*.txt) #check if folder contains txt file to run
	if [ -f "$res" ]; then
                directories[counter]=$f
		let counter=counter+1
	fi
done

for directory in "${directories[@]}"
do	
	for s in $directory/matRad_plan_*.txt
	do

		# Start TOPAS calculation of specific txt file
		echo "Starting calculation for file $(basename "$s") in folder $directory"
		/opt/topas/3.8.1/bin/topas ${s%.*}.txt > ${s%.*}.out > ${s%.*}.err


	done
done

echo "Submitted calculation of ${#directories[@]} folders."

exit 0
EOT








