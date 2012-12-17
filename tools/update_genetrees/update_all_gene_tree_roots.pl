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
use strict;
use warnings;
use JSON;
use DBI;
print "Get DB connection\n";

#my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(-user   => 'treefam_admin',
#                                 -dbname => 'treefam_homology_67hmm',
#                                                                -host   => 'web-mei-treefam',
#                                                                -pass => 'treefam_king1982',
#                                                                -port => '3365');

# MYSQL CONFIG VARIABLES
my $host = "web-mei-treefam";
#my $database = "treefam_homology_67hmm";
my $database = "treefam_production_9_68";
my $user = "treefam_admin";
my $pw = "treefam_king1982";
my $port = "3365";

my $db = DBI->connect("DBI:mysql:database=$database;host=$host;port=$port",$user, $pw) or die "Cannot connect to MySQL server\n";

my $query_gtrt = "select value from gene_tree_root_tag where tag = 'model_name' and root_id = ? ;";
my $query_gtr = "select root_id from gene_tree_root where tree_type ='tree';";

my $update_gtr = "update gene_tree_root set stable_id where root_id = ?;";


my $sth = $db->prepare('select root_id from gene_tree_root where tree_type ="tree"');
my $sql= "select root_id from gene_tree_root where tree_type ='tree'";
my $sql_handle = $db->prepare($sql) or die "Could not prepare\n";
$sql_handle->execute() || die "could not execute\n";

while ( my ($current_id) = $sql_handle->fetchrow_array() ){
	print "\tsearching for $current_id\n";
	#next if $current_id eq '11545517';
	my $sql_gtrt = "select value from gene_tree_root_tag where tag = 'model_name' and root_id = ?";
	my $sql_handle_gtrt = $db->prepare($sql_gtrt) or die "Could not prepare\n";
	$sql_handle_gtrt->execute($current_id) || die "could not execute\n";


	my $model_name;
	while (  ($model_name) = $sql_handle_gtrt->fetchrow_array() ){
		if(!defined($model_name) || $model_name eq '' ){
			die "Could not find model_name for $current_id\n";
		}
		print "Found $model_name\n";
		#print "Update gene_tree_root set stable_id = $model_name where root_id = $current_id\n";
		#exit;
		# prepare template query
		my $sth = $db->prepare("Update gene_tree_root set stable_id = '$model_name' where root_id = $current_id");
		print "Update gene_tree_root set stable_id = '$model_name' where root_id = $current_id\n";
		# execute query with first set of parameters
		my $no_rows = $sth->execute();
		if ($no_rows == 0){
			print  "could not update\n";
		}
		else{
		    print  " updated $no_rows rows\n";
		}
	}

	#last;

}


exit;
my %superhash;

my @families = ("11723527","12255948","12193054","12698788","12506762","11655023","12619130","12405244","11595034");
foreach my $entry (@families){
	
	print "Get Tree object\n";
	$superhash{'treeObject'} = get_tree_object(\%superhash,$db,$entry);
		die "Could not get tree object for $entry\n" if(! exists $superhash{'treeObject'} || !defined($superhash{'treeObject'}));
	print "Got Tree object\n";
	#---------------------------------------------------------------------------------------------------
	print "Get Summary Infos\n";
	if(!get_summary(\%superhash,$superhash{'treeObject'},$entry)){
		die "Could not get sumarrtree object for $entry\n";
	}
	print "Got Summary data\n";
	#---------------------------------------------------------------------------------------------------
	print "Get Sequences\n";
	if(!get_sequences(\%superhash,$superhash{'treeObject'},$entry)){
		die "Could not get sumarrtree object for $entry\n";
	}
	print "Got Sequences\n";
	#---------------------------------------------------------------------------------------------------
	print "Get Tree in Nhx\n";
	if(!get_tree_nhx(\%superhash,$superhash{'treeObject'},$entry)){
		die "Could not get sumarrtree object for $entry\n";
	}
	print "Got Tree in Nhx\n";
	#---------------------------------------------------------------------------------------------------
	print "Get Alignment\n";
	if(!get_alignment(\%superhash,$superhash{'treeObject'},$entry)){
		die "Could not get sumarrtree object for $entry\n";
	}
	print "Got Alignment\n";
	#---------------------------------------------------------------------------------------------------
	print "Get all homologies\n";
	if(!get_homologies(\%superhash,$db,$superhash{'treeObject'},$entry)){
		die "Could not get sumarrtree object for $entry\n";
	}
	print "Got all homologies\n";
	
	
	#print "tree is :".$superhash{'treeNhx'}."\n";
	#print "sequences are\n".$superhash{'sequence_array_json'}."\n";
	#print "homologs\n".$superhash{'homologs_array_json'}."\n";
	
	# sequences
	my $query = "insert into gene_tree_root_tag (root_id,tag,value) values (?, ?, ?) ";
	
	# prepare your statement for connecting to the database
	my $statement = $db->prepare($query);
	# execute your SQL statement
	$statement->execute($entry, 'sequence_array_json', $superhash{'sequence_array_json'});
	
	# alignment
	$statement->execute($entry, 'fasta_aln', $superhash{'fasta_aln'});
	# homologs
	$statement->execute($entry, 'homologs_array_json', $superhash{'homologs_array_json'});
	# tree
	$statement->execute($entry, 'treenhx', $superhash{'treeNhx'});
	
	# no_species
	$statement->execute($entry, 'numspecies', $superhash{'numSpecies'});
}

# update corresponding gene_tree_root entry!

die "Finished analysis\n";


#inserting into gene_tree_root_tag


exit;

sub get_tree_object(){
  my ( $c,$db, $entry ) = (@_);
  my $genetree_adaptor = $db->get_GeneTreeAdaptor;
    my $tree             = $genetree_adaptor->fetch_by_root_id($entry);
    if ( !defined($tree) || $tree eq '' )
    {
        print "Could not get tree\n";
            return;
    }
	return $tree; 
}
sub get_summary(){
  my ( $c,$tree, $entry ) = (@_);

	#get tag-values for this tree
    my $tagvalue_hashref = $tree->get_tagvalue_hash();
    if ( !keys(%$tagvalue_hashref) )
    {
        die "Could not get tagvalue_hashref for tree\n";
    }

	$c->{numSequences}              = $tagvalue_hashref->{gene_count};
    $c->{tree_max_branch}           = $tagvalue_hashref->{tree_max_branch};
    $c->{tree_num_human_peps}       = $tagvalue_hashref->{tree_num_human_peps};
    $c->{tree_num_dup_nodes}        = $tagvalue_hashref->{tree_num_dup_nodes};
    $c->{aln_method}                = $tagvalue_hashref->{aln_method};
    $c->{aln_percent_identity}      = $tagvalue_hashref->{aln_percent_identity};
    $c->{aln_num_residues}          = $tagvalue_hashref->{aln_num_residues};
    $c->{tree_num_spec_nodes}       = $tagvalue_hashref->{tree_num_spec_node};
    $c->{aln_length}                = $tagvalue_hashref->{aln_length};
    $c->{aln_runtime}               = $tagvalue_hashref->{aln_runtime};
    $c->{tree_max_length}           = $tagvalue_hashref->{tree_max_length};
    $c->{buildhmm_runtime_msec}     = $tagvalue_hashref->{buildhmm_runtime_msec};
    $c->{njtree_phyml_runtime_msec} = $tagvalue_hashref->{njtree_phyml_runtime_msec};
    $c->{orthotree_runtime_msec}    = $tagvalue_hashref->{orthotree_runtime_msec};
    $c->{tree_num_leaves}           = $tagvalue_hashref->{tree_num_leaves};

}
sub get_sequences(){
  my ( $c,$tree, $entry ) = (@_);
# Get all sequences
    my $all_leaves = ( $tree->root )->get_all_leaves();
#    $c->log->debug( "Found " . scalar(@$all_leaves) . " members" ) if $c->debug;
#    $c->stash->{region_rows} = $all_leaves;

    my @leaf_array;    
    my %leaf_hash;
    my $counter = 0;
    my @array_of_arrays;
    my $symbol4family = "NaN";
   # get data in array format
    foreach my $leaf(@{$all_leaves}){
        my @tmp_array;
        my $species_name = ($leaf->taxon)->binomial;
        warn "Problem with ".$leaf->taxon->name."\n" if($species_name eq "");
        $species_name =~ s/\s/_/;
        my $species_image_file = "http://localhost:3000/static/images/species_pictures/species_files/thumb_".$species_name.".png";
		push(@tmp_array, "<img src=''$species_image_file'' /> ".($leaf->taxon)->binomial);

        push(@tmp_array,$leaf->stable_id );
        push(@tmp_array,(defined $leaf->display_label ? $leaf->display_label: 'NaN'));
		if(defined $leaf->display_label){ 
		$symbol4family = $leaf->display_label;
			# grep description
			$symbol4family =~ s/-00\d//;
			my $grep_cmd = "grep -w \"$symbol4family\" /Users/fs9/Downloads/mart_export.txt";
			#print "$grep_cmd\n";
			my $grepped_line;# = `$grep_cmd`;
			chomp($grepped_line);
			my ($EnsemblGeneID,$EnsemblProteinID,$HGNCsymbol,$PFAMID,$WikiGeneDescription,$GOTermAccession,$GOTermName,$GOTermDefinition,$GOdomain,$GOTermEvidenceCode,$RfamID) = split("\t",$grepped_line);
			#print "grepped_line: $grepped_line";
			#$c->stash->{treefam}{pfamID} = (defined $PFAMID ? $PFAMID : "No Pfam ID");
			#$c->stash->{treefam}{Wikigenedescription} = (defined $WikiGeneDescription ? $WikiGeneDescription : "No WikiGene description");
			#$c->stash->{treefam}{GOTermAccession} = (defined $GOTermAccession ? $GOTermAccession : "No GOTermAccession");
			#
		}
		push(@tmp_array,$leaf->description );
        #$c->log->debug( "Found line:  ".join(",",@tmp_array)."  members" ) if $c->debug;
        push(@array_of_arrays,\@tmp_array);
    }
    $c->{'sequence_array_json'} =  encode_json \@array_of_arrays;
   
    #$c->log->debug('Family::Tree::get_summary decoded all sequences') if $c->debug; 
    # Count species
    my %species_count;
    #$c->{'numSpecies'} = keys(%species_count);
    foreach(@$all_leaves){ 
         $species_count{($_->taxon)->binomial} += 1;
    }
    $c->{'numSpecies'} = keys(%species_count);
    #$c->log->debug( "Found " . keys(%species_count) . " species" ) if $c->debug;
    return 1;
}
sub get_tree_nhx(){
  my ( $c,$tree, $entry ) = (@_);

  my $root_of_tree = $tree->root;
  $c->{'treeNhx'} = $root_of_tree->nhx_format;
  
  return 1;
  
}
sub get_homologies(){
  my ( $c,$db, $tree, $entry ) = (@_);

#####
    #### Homologies
    ######
    my $HomologyAdaptor = $db->get_HomologyAdaptor;
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
sub get_alignment(){
  my ( $c,$tree, $entry ) = (@_);

# Get Alignment
    my $simpleAlign = $tree->get_SimpleAlign();
    my $fasta_string;
    foreach my $seq ( $simpleAlign->each_seq ){
                $fasta_string .= ">" . $seq->id . "\n" . $seq->seq . "\n";
    }
    $c->{'fasta_aln'} = $fasta_string;
    
  return 1;
  
}
