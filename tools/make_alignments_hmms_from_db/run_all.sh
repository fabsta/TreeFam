for line in $(cat missing_families.txt); 
	do  ./DumpPreviousTreeFamAlnHMMs.pl $line;
done;
