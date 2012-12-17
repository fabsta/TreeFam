#!/bin/bash

NoOfArrayJobs=5
echo About to start $NoOfArrayJobs
NoOfJobsPerArray=250
echo Each job will have $NoOfJobsPerArray
TotalNumberOfJobs=$(cat all_treefamA.txt | wc -l)
echo There are $TotalNumberOfJobs jobs in total

for (( Counter=1; Counter<=TotalNumberOfJobs; Counter=Counter+NoOfJobsPerArray ))
do 
		let lastJob=$(( $Counter + $NoOfJobsPerArray - 1 ))
		echo last Job is $lastJob
		bsub  -o %I.out -e %I.err -J "myArray[$Counter-$lastJob]" /nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/treefam_tools/make_db_quick/make_family_quick.pl -counter input.\$LSB_JOBINDEX
		echo bsub -J "myArray[$Counter-$lastJob]" myJob
done;
