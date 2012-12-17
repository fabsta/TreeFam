#!/bin/bash

NoOfArrayJobs=11
echo About to start $NoOfArrayJobs
NoOfJobsPerArray=500
echo Each job with have $NoOfJobsPerArray
#TotalNumberOfJobs=$(cat first_50_families.txt | wc -l)
TotalNumberOfJobs=5050
echo There are $TotalNumberOfJobs jobs in total

for (( Counter=1; Counter<=TotalNumberOfJobs; Counter=Counter+NoOfJobsPerArray ))
do 
		let lastJob=$(( $Counter + $NoOfJobsPerArray - 1 ))
		echo last Job is $lastJob
		bsub  -o %I.out -e %I.err -J "myArray[$Counter-$lastJob]" /nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/treefam_tools/get_homologs/get_homologs_pairwise.pl -id input.\$LSB_JOBINDEX
		echo bsub -J "myArray[$Counter-$lastJob]" myJob
done;
