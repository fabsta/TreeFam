for line in $(cat helper_files/all_tf9.txt); 
	do  $(pwd)/check_single_copy_genes.pl -id $line;
done;
