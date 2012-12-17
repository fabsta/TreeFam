for line in $(cat all_tf9.txt); 
	do  /nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/treefam_tools/make_db_quick/make_family_quick.pl -id $line;
done;
