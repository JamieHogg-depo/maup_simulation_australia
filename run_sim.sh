#!/bin/bash

for QaS in 1 2 3 4 5 6 7 8 9;
do

	# cur_data
	cur_date=$(date +%Y%m%d)
	nu_zone_seeds=100
	
	# create directories
	mkdir -p CDS_MAUP2/sub_src/$cur_date
	mkdir -p CDS_MAUP2/outputs/$cur_date/lyra_errors
	mkdir -p CDS_MAUP2/outputs/$cur_date/lyra_out
	mkdir -p CDS_MAUP2/outputs/$cur_date/r
	
	# create the unique .sub script files
	specs='QaS'$QaS
	file=$specs'.sub'
	
	# paste the commands in the .sub scripts
	cat>CDS_MAUP2/sub_src/$cur_date/$file<<EOF
#!/bin/bash -l
#PBS -N $specs
#PBS -l ncpus=1
#PBS -l mem=120GB
#PBS -l walltime=80:00:00
#PBS -e CDS_MAUP2/outputs/$cur_date/lyra_errors/$specs
#PBS -o CDS_MAUP2/outputs/$cur_date/lyra_out/$specs

module load r/4.0.3-foss-2020b
module load gdal/3.2.1-foss-2020b

R -e ".libPaths('r_lib'); 
cur_date='$cur_date'; 
QaS=$QaS;
nu_zone_seeds=$nu_zone_seeds;
fit_sa1_data = TRUE;
source('CDS_MAUP2/master.R');"
EOF

	# run each script
		qsub CDS_MAUP2/sub_src/$cur_date/$file
		
done

			