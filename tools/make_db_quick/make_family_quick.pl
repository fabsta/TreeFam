#!/usr/bin/perl
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::Family;
use Bio::EnsEMBL::Compara::DBSQL::GeneTreeAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Compara::GeneTree;
use Bio::SimpleAlign;
use Bio::AlignIO;

use Data::Dumper;
use Scalar::Util;
use Getopt::Long;
use strict;
use warnings;
use JSON;
use File::Temp qw/ tempdir tempfile /;
use TreeFam::Production_modules;

my $entry;
my $counter;
my $switch = "all" ;
my $write_switch = "all" ;
my $write_db = 1;
my $clean_old_entries  = 1;
my $db_host = "web-mei-treefam";
my $db_database = "treefam_production_9_68";
#my $db_database = "treefam_homology_67hmm";
my $db_user = $ENV{'TREEFAM_DB_USER'};
my $db_pass = $ENV{'TREEFAM_DB_PASS'};
my $db_port = 3365;
my $result = GetOptions ("id=s" => \$entry,    # numeric
                      "switch=s"   => \$switch,      # string
                      "write_switch=s"  => \$write_switch,
                        "counter=s" => \$counter);

if(!defined($entry) && $counter){
    print "checking counter $counter\n";
    ## get entry 
    if($counter =~ /input/){
        $counter =~ s/input\.//;
    }
    my @all_jobs = `cat all_tf9.txt`; 
    #print "will try $treefamFamilyID_number\n";
    $entry = $all_jobs[$counter-1];
    chomp($entry);
    }
die "no valid entry $entry supplied" if !($entry =~ /\d+/); 
print "parameters are: $switch, $write_switch, $entry\n";
#exit;

my $mappings_dir = "/lustre/scratch109/sanger/fs9/";

our $description_file = "$mappings_dir/mappings/ensembl_biomart/biomart_descriptions.txt";
our $go_file = "$mappings_dir/mappings/ensembl_biomart/biomart_go.txt";
our $hgnc_file = "$mappings_dir/mappings/hgnc/hngc_mappings.txt";
our $uniprot_file = "$mappings_dir/mappings/uniprot/ens2all_mapping.txt";
our $pfam_file = "$mappings_dir/mappings/pfam/pfam_predictions.txt";
our $previous_treefam_file = "$mappings_dir/mappings/previous_treefam/familyA.txt.table";
our $wikigene_file = "$mappings_dir/mappings/wikigene/results_wikigene_hsapiens_gene_ensembl";
our $hmmer_scores_file = "$mappings_dir/mappings/hmmer_scores/hmmer_scores.txt";
our $treebest_binary = "/nfs/users/nfs_f/fs9/bin/treebest-1.9.2_i686-linux/treebest";

print "Check existence of mapping files";
die "File desc missing\n" if(!-e $description_file|| !-s $description_file );
die "File go missing\n" if(!-e $go_file|| !-s $go_file );
die "File hgnc missing\n" if(!-e $hgnc_file|| !-s $hgnc_file );
die "File pfam missing\n" if(!-e $pfam_file|| !-s $pfam_file );
die "File hmmer-scores missing\n" if(!-e $hmmer_scores_file || !-s $hmmer_scores_file);
die "File uniprot missing\n" if(!-e $uniprot_file || !-s $uniprot_file);

print "Get DB connection\n";


### Connect to DB
my $registry = 'Bio::EnsEMBL::Registry';
my $registry_file = '/nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/production_treefam_reg_conf.pl';
#$ENV{'ENSEMBL_REGISTRY'} = '/nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/production_treefam_reg_conf.pl';
TreeFam::Production_modules::get_db_connection({"registry" => $registry,"registry_file" => $registry_file});
die "Could not connect to DB. Check registry in $registry_file\n" if !$registry;


my $member_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'Member' );
if(!$member_adaptor){
    die "Could not load member_adaptor\n";
}
#my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(-user   => '$db_user',
                                                    #-dbname => '$db_database',
                                                    #-host   => 'web-mei-treefam',
                                                                #-pass => '$db_pass',
                                                                #-port => '$db_port');
my $db;
my %superhash;
my %ext_counts;
#my @families = ("12619130");
#my @families = ("12698788","12391116");  # TF101005
##my @families = ("12391116"); # TF101024
#my @families = ("12255948"); # TF101003
#my @families = ("12506762"); # TF101006
#my @families = ("11975910"); # TF101001
#my @families = ("12663583"); # TF101001
#my @families = ("11922233"); # TF101001

#my @all_jobs = `cat treefam8_allfamilies3.txt`; 
#my $treefamFamilyID_number       = $ARGV[0];
#$treefamFamilyID_number =~ s/input\.//;
#print "will try $treefamFamilyID_number\n";
#my $treefamFamilyID = $all_jobs[$treefamFamilyID_number-1];
#chomp($treefamFamilyID);
#print "will try $treefamFamilyID\n";

# causing problem: "11975910" "12255948"  

#my @families = ("12619130","12255948", "11975910","12698788","12506762","11655023","12405244","11595034");
# todo: replace with iterator over all families
#TREE:
#foreach my $entry (@families){
    my %sequences_hash;	
    my %species_count;	
    my $mark_family = 0;
    die "wait, entry is invalid: $entry\n" if !defined($entry) || $entry eq '';
	print "Looking at: $entry\n";
    print "\tGet Tree object\n";

	$superhash{'treeObject'} = TreeFam::Production_modules::get_tree_object({ "registry" => $registry,"entry" =>$entry});

    if(! exists $superhash{'treeObject'} || !defined($superhash{'treeObject'})){
	    warn "\tCould not get tree object for $entry\n"; 
        exit;#next TREE;
    }
    my ($treefam_id, $treefam_name) = ($entry,$superhash{'treeObject'}->stable_id);
	print "\tGot Tree object for $entry\n";
        
    if($clean_old_entries){
            TreeFam::Production_modules::clean_old_entries({
                    "entry" => $entry,
                    "db_host" => $db_host,
                    "db_user" => $db_user,
                    "db_pass" => $db_pass,
                    "db_port" => $db_port,
                    "db_database" => $db_database
            });
    }
    #---------------------------------------------------------------------------------------------------
# SUMMARY
#---------------------------------------------------------------------------------------------------

    print "Get Summary Infos\n";
	if(!TreeFam::Production_modules::get_summary({
                                    "superhash" => \%superhash, 
                                    "treeObject" => $superhash{'treeObject'},
                                    "entry" => $entry})){
		                                die "Could not get summary object for $entry\n";
	}
    if(!exists $superhash{tree_num_leaves} || $superhash{tree_num_leaves} == 0){
        warn "No number of sequences, but we can set this value later \n";
        $mark_family = 1;
    }
	print "Got Summary data\n";
#---------------------------------------------------------------------------------------------------
# SEQUENCES
#---------------------------------------------------------------------------------------------------
if($switch eq "all" || $switch eq "sequences" || $switch eq "tree"){	
	print "Get Sequences\n";
	if(!TreeFam::Production_modules::get_sequences({
                                    "superhash" => \%superhash,
                                    "sequences_hash" => \%sequences_hash,
                                    "treeObject" => $superhash{'treeObject'},
                                    "species_count" =>  \%species_count})){
		                                    die "Could not get sequences object for $entry\n";
	}
    print "Got Sequences\n";
    ## In case we re-use a tree
    # we have to set the number of sequences
    if(!exists $superhash{tree_num_leaves}){
        $superhash{tree_num_leaves} = keys(%sequences_hash);
        print "TREE_NUM_LEAVES set to ".keys(%sequences_hash)."\n"; 
    }
}
#---------------------------------------------------------------------------------------------------
# ANNOTATE SEQUENCES
#---------------------------------------------------------------------------------------------------
if($switch eq "all" || $switch eq "sequences" || $switch eq "tree" || $switch eq "family"){	
	print "Get sequence annotation\n";
    if(!TreeFam::Production_modules::get_sequence_annotations({
                                    "sequences_hash" => \%sequences_hash,
                                    "sequence_array_json" => \$superhash{'sequence_array_json'},
                                    "hgnc_file" => $hgnc_file,
                                    "wikigene_file" => $wikigene_file,
                                    "pfam_file" => $pfam_file,
                                    "hmmer_scores_file" => $hmmer_scores_file,
                                    "uniprot_file" => $uniprot_file,
                                    "ext_counts" => \%ext_counts,
                                    "db_adaptor" => $member_adaptor
                                    })){
		                                    die "Could not get annotation object for $entry\n";
	}
	print "Got all annotations\n";
	#print "sequences are\n".$superhash{'sequence_array_json'}."\n";
}

#---------------------------------------------------------------------------------------------------
# FAMILY ANNOTATION 
#---------------------------------------------------------------------------------------------------
if($switch eq "all" || $switch eq "family" ){	
	print "Get family annotation\n";
    if(!TreeFam::Production_modules::get_family_annotations({
                                        "superhash" => \%superhash,
                                        "sequences_hash" => \%sequences_hash,
                                        "treeObject" =>  $superhash{'treeObject'}, 
                                        "species_count" => \%species_count,
                                        "previous_treefam_file" => $previous_treefam_file})){
		                                    die "Could not get family annotation object for $entry\n";
	}
	print "Got all annotations\n";
	#print "sequences are\n".$superhash{'sequence_array_json'}."\n";
}

#---------------------------------------------------------------------------------------------------
# TREE (NEWICK)
#---------------------------------------------------------------------------------------------------
if($switch eq "all" || $switch eq "tree" || $switch eq "onlyTree"){	
	print "Get Tree in Nhx\n";
	if(!TreeFam::Production_modules::get_tree_nhx({
                                        "superhash" => \%superhash,
                                        "treeObject" => $superhash{'treeObject'},
                                        "entry" =>  $entry })){
	                                            	die "Could not get newick tree object for $entry\n";
	}
	print "Got Tree in Nhx\n";
    #print "here it is: ".$superhash{'treeNhx'}."\n";
    #next TREE;
#---------------------------------------------------------------------------------------------------
# Make tree image
#---------------------------------------------------------------------------------------------------
	print "Get Tree image \n";
	#if(!get_tree_image(\%superhash,$superhash{'treeNhx'}, $superhash{"tree_num_leaves"})){
	#	warn "Could not get image for $entry\n";
	#}
	print "Got Tree image \n";

#---------------------------------------------------------------------------------------------------
# TREE (Phyloxml)
#---------------------------------------------------------------------------------------------------
	print "Get Tree in phyloxml\n";
    if(!TreeFam::Production_modules::get_tree_phyloxml({
                                                "superhash" => \%superhash,
                                                "registry" => $registry,
                                                "treeNhx" => $superhash{'treeNhx'}, 
                                                "sequences_hash" => \%sequences_hash, 
                                                "entry" => $entry, 
                                                "treefam_name"=> $treefam_name})){
		                                                warn "Could not get phyloxml object for $entry\n";
	}
	print "Got Tree in phyloxml\n";
    print $superhash{'phyloxml'};
    #die "stop here for now\n";
}
#---------------------------------------------------------------------------------------------------
# ALIGNMENT
#---------------------------------------------------------------------------------------------------
if($switch eq "all" || $switch eq "alignment"){	
	print "Get Alignment\n";
    unless($mark_family){
        if(!TreeFam::Production_modules::get_alignment({
                                        "superhash" => \%superhash,
                                        "treeObject" => $superhash{'treeObject'}, 
                                        "entry" => $entry})){
	                                                warn "Could not get alignment object for $entry\n";
	    }
    }
	print "Got Alignment\n";
}
#---------------------------------------------------------------------------------------------------
# HOMOLOGIES
#---------------------------------------------------------------------------------------------------
if($switch eq "all" || $switch eq "homologies"){	
	print "Get all homologies\n";
	#if(!get_homologies(\%superhash,$registry,$superhash{'treeObject'},$entry)){
	#	die "Could not get sumarrtree object for $entry\n";
#}
	print "Got all homologies\n";
}

#---------------------------------------------------------------------------------------------------
# Checking if all information is correct
#---------------------------------------------------------------------------------------------------


	#print "tree is :".$superhash{'treeNhx'}."\n";
	#print "sequences are\n".$superhash{'sequence_array_json'}."\n";
	#print "homologs\n".$superhash{'homologs_array_json'}."\n";
    #exit;	
	# sequences
#---------------------------------------------------------------------------------------------------
#  Writing information to db 
#---------------------------------------------------------------------------------------------------
if($write_db){

    my $query = "insert into gene_tree_root_tag (root_id,tag,value) values (?, ?, ?) ";
	
	# prepare your statement for connecting to the database
    my $genetree_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'GeneTree' );
    if(!$genetree_adaptor){
        die "Could not load genetree_adaptor\n";
    }
    my $statement = $genetree_adaptor->prepare($query);
	
    # SEQUENCES
    if($write_switch eq "all" || $write_switch eq "write_sequences"){
        print "Writing sequence_array_json\n";
        $statement->execute($entry, 'sequence_array_json', $superhash{'sequence_array_json'});
    }
    # ANNOTATIONS 
    if($write_switch eq "all" || $write_switch eq "write_family"){
        # Fam_symbol
        if(!$superhash{'fam_symbol'}){warn "DB_ERROR: Could not save fam_symbol\n";}
        else{$statement->execute($entry, 'fam_symbol', $superhash{'fam_symbol'});}
        # Fam_description
        if(!$superhash{'fam_description'}){warn "DB_ERROR: Could not save fam_symbol\n";}
        else{$statement->execute($entry, 'fam_description', $superhash{'fam_description'});}
        # Fam_n_seed
        if(!$superhash{'fam_n_seed'}){warn "DB_ERROR: Could not save fam_n_seed\n";}
        else{$statement->execute($entry, 'fam_n_seed', $superhash{'fam_n_seed'});}
        # Fam_n_seed
        if(!$superhash{'fam_n_seed'}){warn "DB_ERROR: Could not save fam_n_seed\n";}
        else{$statement->execute($entry, 'fam_n_full', $superhash{'fam_n_full'});}
        # NumSequences
        if(!$superhash{'numSequences'}){warn "DB_ERROR: Could not save numSequences\n";}
        else{$statement->execute($entry, 'numSequences', $superhash{'numSequences'});}
        # hngc
        if(!$superhash{'hgnc'}){warn "DB_ERROR: Could not save hgnc\n";}
        else{$statement->execute($entry, 'hgnc', $superhash{'hgnc'});}
        # pfam
        if(!$superhash{'pfam'}){warn "DB_ERROR: Could not save pfam\n";}
        else{$statement->execute($entry, 'pfam', $superhash{'pfam'});}
        # numSpecies
        if(!$superhash{'numSpecies'}){warn "DB_ERROR: Could not save numSpecies\n";}
        else{$statement->execute($entry, 'numSpecies', $superhash{'numSpecies'});}
        # wikigene
        if(!$superhash{'wikigene'}){warn "DB_ERROR: Could not save wikigene\n";}
        else{$statement->execute($entry, 'wikigene', $superhash{'wikigene'});}
        # taxa_count
        if(!$superhash{'taxa_count'}){warn "DB_ERROR: Could not save taxa_count\n";}
        else{$statement->execute($entry, 'taxa_count', $superhash{'taxa_count'});}
    } 

    
    # alignment
    if($write_switch eq "all" || $write_switch eq "write_alignment"){
        if(!$superhash{'fasta_aln'}){
            warn "DB_ERROR: Could not save alignment\n";
        }
        else{
            $statement->execute($entry, 'fasta_aln', $superhash{'fasta_aln'});
        }
    }
    # homologs
    #if($write_switch eq "all" || $write_switch eq "write_homologs"){
	#    $statement->execute($entry, 'homologs_array_json', $superhash{'homologs_array_json'});
    #}
    # tree
    if($write_switch eq "all" || $write_switch eq "write_tree"){
	    if(!$superhash{'treeNhx'}){
            warn "DB_ERROR: Could not save tree (NHX)\n";
        }
        else{
            $statement->execute($entry, 'treenhx', $superhash{'treeNhx'});
        }
    }

    # tree image
    if($write_switch eq "all" || $write_switch eq "write_tree_image"){
        if(!$superhash{'tree_image_png'}){
            warn "DB_ERROR: Could not save tree image\n";
        }
        else{
            $statement->execute($entry, 'tree_image_png', $superhash{'tree_image_png'});
        }
    }
    # tree phyloxml
    if($write_switch eq "all" || $write_switch eq "write_tree_phyloxml"){
        if(!$superhash{'treephyloxml'}){
            warn "DB_ERROR: Could not save tree (phyloxml)\n";
        }
        else{
            $statement->execute($entry, 'treephyloxml', $superhash{'treephyloxml'});
        }
    }
    # write_sequence 
    if($write_switch eq "all" || $write_switch eq "write_sequence_mappings"){
        if(!keys(%sequences_hash)){
            warn "DB_ERROR: Could not save sequence mappings\n";
        }
        else{
            #print Dumper %ext_counts;
            print "well, this was before calling the function\n";
            &write_sequence_mappings({ "sequence_href" => \%sequences_hash, "treefam_tree_id" => $treefam_id,"treefam_tree_name" => $treefam_name, "ext_counts" => \%ext_counts});
        }
    }

# no_species
    }
#die "looked at last one (maybe 1) family\n";
#}

# update corresponding gene_tree_root entry!

die "Finished analysis\n";


#inserting into gene_tree_root_tag


exit;


sub get_tree_image(){
    my ( $c,$newickTree, $no_sequences ) = (@_);
    #print Dumper $newickTree;
    #print Dumper $no_sequences;
    my ($fh, $tree_file) = tempfile( );
    if(!$no_sequences ){
        warn "Could not determine  no of sequences. Skipping making an image\n";
        return 0;
    }
    print "Tree to file $tree_file with $no_sequences\n";
    
    # save newick string as file
    open my $nw_tree_out, ">", $tree_file or die "Could not open $tree_file\n";
    print {$nw_tree_out} $newickTree."\n";
    close $nw_tree_out || die "Could not close $tree_file\n";
    if(!-e $tree_file || ! -s $tree_file){
        warn "Problem saving tree to file $tree_file\n";
        return 0;
    } 
    my $outtree = $tree_file.".eps";
    my $tree_png = $tree_file.".png";

    # read the number of sequences in file $mfa9:
    # calculate the height of the image:
    my $height = $no_sequences * 11;
    print "$treebest_binary export $tree_file -y $height > $outtree";
    #exit;
    system "$treebest_binary export $tree_file -y $height > $outtree";
    print "convert $outtree $tree_png\n";	
    system "convert $outtree $tree_png";	
    
    $c->{"tree_image_png"} = `cat $tree_png` ;
    #exit;
    return -e $tree_png; 

}


sub get_homologies(){
  my ( $c,$registry, $tree, $entry ) = (@_);

#####
    #### Homologies
    ######
    my $HomologyAdaptor = $registry->get_HomologyAdaptor;
    my $homologies      = $HomologyAdaptor->fetch_all_by_tree_node_id( ($tree->root)->node_id );
    print "Found " . scalar(@$homologies) . " homologies \n";
    my %homology_type_hash;
    my @per_species_homologs;
    my %taxonomy_count;
    my $totalHomologs = 0;
    my @array_of_homologs_array;
    foreach my $this_homology (@$homologies)
    {
        $homology_type_hash{ $this_homology->description }++;
        my $taxonomy_level = $this_homology->taxonomy_level();
        $taxonomy_count{$taxonomy_level}{ $this_homology->description }++;
        #next;
        my $all_members = $this_homology->get_all_Members();
        # Assuming that each homology consists of two members only
        my ($member1,$member2) = ($all_members->[0],$all_members->[1]); 
        $per_species_homologs[$totalHomologs]{'type'} =  $this_homology->description;
        $per_species_homologs[$totalHomologs]{'member1'} =  $member1;
        $per_species_homologs[$totalHomologs]{'member2'} =  $member2;
        my $species1_name = ($member1->taxon)->binomial;
        $species1_name =~ s/\s/_/;
        my $species2_name = ($member2->taxon)->binomial;
        $species2_name =~ s/\s/_/;
        
        my $species_image_file = "http://localhost:3000/static/images/species_pictures/species_files/thumb_".$species1_name.".png";
        my $species2_image_file = "http://localhost:3000/static/images/species_pictures/species_files/thumb_".$species2_name.".png";
        my @tmp_array;
        push(@tmp_array, "<img src=''$species_image_file'' /> ".($member1->taxon)->binomial);
        push(@tmp_array,$member1->stable_id );
        my $translated_type;
        if($this_homology->description eq 'ortholog_one2one'){ $translated_type = "1-1 Orthologs";}
        if($this_homology->description eq 'within_species_paralog'){ $translated_type = "Inparalogs";}
        if($this_homology->description eq 'ortholog_one2many'){ $translated_type = "1-many orthologs";}
        if($this_homology->description eq 'possible_ortholog'){ $translated_type = "possible orthologs";}
        if($this_homology->description eq 'apparent_ortholog_one2one'){ $translated_type = "apparent 1-1 Orthologs";}
        #if($translated_type eq ''){$c->log->debug("no type for   ".$this_homology->description." ") if $c->debug;}
        push(@tmp_array, $translated_type);
        push(@tmp_array, "<img src=''$species2_image_file'' /> ".($member2->taxon)->binomial);
        push(@tmp_array,$member2->stable_id);
        push(@array_of_homologs_array,\@tmp_array); 
        $totalHomologs++; 
   }
   print "finished reading homologies\n";
   $c->{"homologs_array_json"} = encode_json \@array_of_homologs_array; 
	return 1;
}






sub read_description_file {
my ($file,$prot_hashref,) = (@_);
    open my $OPEN_FH,"<",$file or die "Could not open $file\n";
    my $first_line = <$OPEN_FH>;
    foreach my $line (<$OPEN_FH>){
        my ($counter,$gene_id,$protein_id,$description,$ensembl_family) = split(" ",$line);
        $protein_id =~ s/"//g;
        $description =~ s/"//g;
        $ensembl_family =~ s/"//g;
        
        print "$counter,$gene_id,$protein_id,$description,$ensembl_family \n";
        
        $prot_hashref->{$protein_id}{'description'} = $description;
        $prot_hashref->{$protein_id}{'ensembl_family'} = $ensembl_family;
        last;
    }
    close $OPEN_FH  or die "Could not close $file\n";

}


sub get_description_hits {
    my ($id) = (@_);
    my $file_to_search = $description_file;
    my $grep_command = "grep $id $file_to_search";
    my @result_lines = `$grep_command`;

### Parse
    
    return ;
}


sub read_go_file {
    my ($file,$prot_hashref,) = (@_);
    open my $OPEN_FH,"<",$file or die "Could not open $file\n";
    my $first_line = <$OPEN_FH>;
    foreach my $line (<$OPEN_FH>){
        my ($counter, $protein_id,$go_id,$go_name, $go_definition, $go_evidence, $go_namespace) = split(" ",$line);
        $protein_id =~ s/"//g;
        $go_id =~ s/"//g;
        $go_name =~ s/"//g;
        $go_definition =~ s/"//g;
        $go_evidence =~ s/"//g;
        $go_namespace =~ s/"//g;
       print "$counter, $protein_id,$go_id,$go_name, $go_definition, $go_evidence, $go_namespace\n"; 
        $prot_hashref->{$protein_id}{'go_id'} = $go_id;
        $prot_hashref->{$protein_id}{'go_name'} = $go_name;
        $prot_hashref->{$protein_id}{'go_definition'} = $go_definition;
        $prot_hashref->{$protein_id}{'go_evidence'} = $go_evidence;
        last;
        }
    close $OPEN_FH  or die "Could not close $file\n";

}



 
# now write the output
 
exit;
#die "finished looking at nodes\n";



sub write_sequence_mappings{
	my ($arg_ref) = @_;
	my $sequence_href = $arg_ref->{sequence_href};
	my $treefam_tree_id = $arg_ref->{treefam_tree_id};
	my $treefam_tree_name = $arg_ref->{treefam_tree_name};
	my $ext_counts = $arg_ref->{ext_counts};
    #print Dumper $ext_counts;
    my %family_level_hits;
    my $dbh = DBI->connect("dbi:mysql:$db_database;host=$db_host:$db_port","$db_user","$db_pass") or die "Connection Error: $DBI::errstr\n";
    my $extID2seq_insert = 'insert into xrefID2Sequence (external_db_id,external_db_id_name,db,member_id,gene_tree_id,gene_tree_stable_id) values (?,?,?,?,?,?)';
    my $insert_handle = $dbh->prepare_cached('insert into xrefID2Sequence (external_db_id,external_db_id_name,db,member_id,gene_tree_id,gene_tree_stable_id, description) values (?,?,?,?,?,?,?)');
    SEQ:
    foreach my $seq(keys(%{$sequence_href})){
           # ge
           my $member_id = $sequence_href->{$seq}{'member_id'};
           if($member_id eq ''){
                    warn "Could not save mappings for $seq\n"; 
                    warn Dumper $seq;
                    next SEQ;
                    }
    # PFAM
        if(keys(%{$sequence_href->{$seq}{'pfam_hits'}})){
            PFAMHIT:
            foreach my $pfam_hit(keys(%{$sequence_href->{$seq}{pfam_hits}})){
            my $ratio = $ext_counts->{"pfam_counts"}{$pfam_hit} / keys(%{$sequence_href});
                my ($id,$name, $a_start,$a_end,$evalue)  =  ($sequence_href->{$seq}{pfam_hits}{$pfam_hit}{hmm_name},
                                                            $sequence_href->{$seq}{pfam_hits}{$pfam_hit}{description},
                                                            $sequence_href->{$seq}{pfam_hits}{$pfam_hit}{alignment_start},
                                                            $sequence_href->{$seq}{pfam_hits}{$pfam_hit}{alignment_end},
                                                            $sequence_href->{$seq}{pfam_hits}{$pfam_hit}{evalue} );
            if($ratio < 0.05){
                print "skipped $id -> $pfam_hit ($ratio ".$ext_counts->{"pfam_counts"}{$pfam_hit}." / ".keys(%{$sequence_href}).")\n";
                next PFAMHIT;
            }
                # lets get rid of pfam version numbers
                $id =~ /(PF\d+)/;
                $id = $1;
                my $description = "From: $a_start To: $a_end evalue: $evalue";
               $insert_handle->execute($id,$name,"Pfam",$member_id,$treefam_tree_id,$treefam_tree_name,$description);
               #or warn "Could not insert: $id,$name,\"Pfam\",$member_id,$treefam_tree_id,$treefam_tree_name\n";
                $family_level_hits{'pfam'}{$name} = $id;
                }
        }
    # UniProt
        if(keys(%{$sequence_href->{$seq}{'uniprot_hits'}})){
            foreach my $uniprot_hit_db(keys(%{$sequence_href->{$seq}{uniprot_hits}})){
                my $ext_id = $sequence_href->{$seq}{uniprot_hits}{$uniprot_hit_db}; 
                $insert_handle->execute($ext_id,"NULL",$uniprot_hit_db,$member_id,$treefam_tree_id,$treefam_tree_name,"");
                #or warn "Could not insert: $ext_id,\"NULL\",$uniprot_hit_db,$member_id,$treefam_tree_id,$treefam_tree_name\n";
            }
        }
        # HGNC
        if(keys(%{$sequence_href->{$seq}{'hgnc_hits'}})){
            foreach my $hgnc_name(keys(%{$sequence_href->{$seq}{hgnc_hits}})){
                my $hgnc_id = $sequence_href->{$seq}{hgnc_hits}{$hgnc_name};
                if($hgnc_name eq "NaN" || !defined($hgnc_name) ){next;} 
                $insert_handle->execute($hgnc_id,$hgnc_name,"HGNC",$member_id,$treefam_tree_id,$treefam_tree_name,"");
                #or warn "Could not insert: $hgnc_id,$hgnc_name,\"HGNC\",$member_id,$treefam_tree_id,$treefam_tree_name\n";
                $family_level_hits{'hgnc'}{$hgnc_name} = $hgnc_id;
            }
        }
        if(keys(%{$sequence_href->{$seq}{'wikigene_hits'}})){
            foreach my $wikigene_hit(keys(%{$sequence_href->{$seq}->{wikigene_hits}})){
                my ($id,$description)  = ($sequence_href->{$seq}->{wikigene_hits}{$wikigene_hit}{id},$sequence_href->{$seq}->{wikigene_hits}{$wikigene_hit}{description});
                $insert_handle->execute($wikigene_hit,$id,"Wikigene",$member_id,$treefam_tree_id,$treefam_tree_name,$description);
                #or warn "Could not insert: $wikigene_hit,\"\",\"WIKIGENE\",$member_id,$treefam_tree_id,$treefam_tree_name\n";
                $family_level_hits{'wikigene'}{$wikigene_hit}{'id'} = $id;
                $family_level_hits{'wikigene'}{$wikigene_hit}{'description'} = $description;
            }
        }
    }
    print "Here are the summaries for the family table\n";
    # to save: hgnc, pfam
     #  BRCA2          | 675                 | HGNC |     11990942 | TF105041            | NULL
    my $insert_fam_handle = $dbh->prepare_cached('insert into xrefID2Family (external_db_id,external_db_id_name,db,gene_tree_id,gene_tree_stable_id, description) values (?,?,?,?,?,?)');
    foreach my $db(keys(%family_level_hits)){
        foreach my $key(keys(%{$family_level_hits{$db}})){
                my ($value) = $family_level_hits{$db}{$key};
                my $insert_string;
                if($db eq 'wikigene'){
                    my ($id,$description)  = ($family_level_hits{$db}{$key}{id}, $family_level_hits{$db}{$key}{description});
                    #$insert_string = "$key, $id,$db,$treefam_tree_id,$treefam_tree_name,$description";
                    print "inserting: $key, $id,$db,$treefam_tree_id,$treefam_tree_name,\"$description\"\n";
                    $insert_fam_handle->execute($key, $id,$db,$treefam_tree_id,$treefam_tree_name,$description);
                }
                else{
                    #$insert_string = "$key, $value,$db,$treefam_tree_id,$treefam_tree_name,\"\"";
                    print "inserting: $key, $value,$db,$treefam_tree_id,$treefam_tree_name,\"\"\n";
                    $insert_fam_handle->execute($key, $value,$db,$treefam_tree_id,$treefam_tree_name,"");
                    }
        }
    }
}

