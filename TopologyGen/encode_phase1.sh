#! /bin/bash
for (( kk=100; kk<= 100; kk=kk+1)) #for the test version
do
	for (( ii=0; ii<=29; ii=ii+1)) #for the camera index	
	do	
		echo "
		### Job name
		#PBS -N Encode"$ii"

		### out files
		#PBS -e ./log/Encode"$ii".err
		#PBS -o ./log/Encode"$ii".log

		### put the job to which queue (qwork)
		#PBS -q qwork
		" > ./tempFile/Encode$ii.sh

		echo '
		# show the time and information
		echo Working directory is $PBS_O_WORKDIR
		cd $PBS_O_WORKDIR
		echo Running on host `hostname`
		echo Start time is `date`
		time1=`date +%s`
		echo Directory is `pwd`
		' >> ./tempFile/Encode$ii.sh

		echo "
		# run the script
		MCR_tmp=\"/tmp/mcr_\$RANDOM\"
		mkdir \$MCR_tmp
		export MCR_CACHE_ROOT=\$MCR_tmp
		python encode.py --phase 1 --camera "$ii"
		rm -rf \$MCR_tmp                                                
		" >>  ./tempFile/Encode$ii.sh

		echo '
		echo End time is `date`
		time2=`date +%s`
		echo Computing time is `echo $time2-$time1 | bc` sec
		' >>   ./tempFile/Encode$ii.sh
	
		qsub ./tempFile/Encode$ii.sh
	done
done
