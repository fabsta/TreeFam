#!/bin/bash

NoOfArrayJobs=100
echo About to start $NoOfArrayJobs
NoOfJobsPerArray=150
echo Each job with have $NoOfJobsPerArray
TotalNumberOfJobs=$(cat treefam8_allfamilies3.txt | wc -l)
echo There are $TotalNumberOfJobs jobs in total

for (( Counter=1; Counter<=TotalNumberOfJobs; Counter=Counter+NoOfJobsPerArray ))
do 
		let lastJob=$(( $Counter + $NoOfJobsPerArray - 1 ))
		echo last Job is $lastJob
		bsub  -o %I.out -e %I.err -J "myArray[$Counter-$lastJob]" /nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/treefam_tools/get_HMMALING/myjob.pl input.\$LSB_JOBINDEX
		echo bsub -J "myArray[$Counter-$lastJob]" myJob
done;
