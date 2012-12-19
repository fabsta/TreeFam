#created
#===============================================================================
#
#         FILE: HomologyHelper.pm
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Fabian Schreiber (), fs9@sanger.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 09/26/2012 16:15:19
#     REVISION: ---
#===============================================================================
package TreeFam::Production_modules;

use strict;
use warnings;
use Bio::Phylo::Factory;
use Bio::Phylo::IO qw'parse unparse';
use Bio::Phylo::Util::CONSTANT qw':objecttypes :namespaces';
use Bio::TreeIO;
use Bio::Phylo::Forest::Tree;
use Data::Dumper;
use Bio::EnsEMBL::Registry; 
use File::Temp qw/ tempdir tempfile /;
use JSON;
sub get_db_connection{
	my ($arg_ref) = @_;
	my $registry = $arg_ref->{registry};
	my $registry_file = $arg_ref->{registry_file};
    
    my $reg = Bio::EnsEMBL::Registry->load_all($registry_file);
    return defined($reg)? $reg : undef;
}
sub get_tree_object{
	my ($arg_ref) = @_;
	my $registry = $arg_ref->{registry};
	my $entry = $arg_ref->{entry};
    
    my $genetree_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'GeneTree' );
    if(!$genetree_adaptor){
        die "Could not load genetree_adaptor\n";
    }

    #my $genetree_adaptor = $registry->get_GeneTreeAdaptor;
    my $tree             = $genetree_adaptor->fetch_by_root_id($entry);
    if ( !defined($tree) || $tree eq '' )
    {
        print "Could not get tree\n";
            return;
    }
	return $tree; 
}

sub clean_old_entries{
	my ($arg_ref) = @_;
	my $entry = $arg_ref->{entry};
	my $db_host = $arg_ref->{db_host};
	my $db_user = $arg_ref->{db_user};
	my $db_pass = $arg_ref->{db_pass};
	my $db_port = $arg_ref->{db_port};
	my $db_database = $arg_ref->{db_database};
 
### do we need to clean old entry?
            #my $treeObject = $superhash{'treeObject'};
            my @tags_to_delete = ("fam_description", "fam_n_full", "hgnc", "fam_n_seed", "fam_full_description","fam_symbol", "fasta_aln", "homologs_array_json","numspecies","pfam","treenhx","treephyloxml","tree_image_png", "sequence_array_json","taxa_count","wikigene", "numSequences");
            #foreach my $tag(@tags_to_delete){
                #print "Deleting tag $tag ...";
                ##print "\n";
                ##exit;
                #my $remove_success = $treeObject->remove_tag($tag);
                #print "".(($remove_success)? "yes":"no")."\n"; 
            #}
            foreach my $tag(@tags_to_delete){
                my $mysql_command = "mysql -h $db_host -P $db_port -u$db_user -p$db_pass -e 'DELETE FROM $db_database.gene_tree_root_tag where root_id = $entry and tag =\"$tag\";'";
                print $mysql_command."\n";
                system($mysql_command);
            }
            # Now delete entries from xrefID2Sequence, 
            my $mysql_command_xrefID2Sequence = "mysql -h $db_host -P $db_port -u$db_user -p$db_pass -e 'DELETE FROM $db_database.xrefID2Sequence where gene_tree_id =  $entry;'";
                print $mysql_command_xrefID2Sequence."\n";
                system($mysql_command_xrefID2Sequence);
            # Now delete entries from xrefID2Sequence, 
            my $mysql_command_xrefID2family = "mysql -h $db_host -P $db_port -u$db_user -p$db_pass -e 'DELETE FROM $db_database.xrefID2Family where gene_tree_id =  $entry;'";
            print $mysql_command_xrefID2family."\n";
            system($mysql_command_xrefID2family);
            
}
sub get_summary{
	my ($arg_ref) = @_;
	my $tree = $arg_ref->{treeObject};
	my $c = $arg_ref->{superhash};
	my $entry = $arg_ref->{entry};

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
    return 1;
}
sub get_sequences{
	my ($arg_ref) = @_;
	my $c = $arg_ref->{superhash};
	my $sequence_hash = $arg_ref->{sequences_hash};
	my $tree = $arg_ref->{treeObject};
	my $species_count = $arg_ref->{species_count};
	
# Get all sequences
    my $all_leaves = ( $tree->root )->get_all_leaves();
    my @leaf_array;    
    my %leaf_hash;
    my $counter = 0;
    my @array_of_arrays;
    my $symbol4family = "NaN";
   # get data in array format
    foreach my $leaf(@{$all_leaves}){
        my @tmp_array;
        my $prot_seq_id = $leaf->stable_id;
        ### Save each sequence with 
        # use member id to search hmmer_scores 
        $sequence_hash->{$prot_seq_id}{"member_id"} = $leaf->dbID;
        # gene member id
        my $gene4prot = $leaf->gene_member;
        $sequence_hash->{$prot_seq_id}{"geneID"} = $gene4prot->stable_id;
        # sequence object
        #$sequence_hash->{$prot_seq_id}{"sequence_object"} = $leaf;
        my $species_name = ($leaf->taxon)->binomial;
        if(!defined($species_name) || $species_name eq ""){
            warn "Problem with ".$leaf->taxon->name."\n";
            $species_name = $leaf->taxon->name;
        }
        # taxon ID
        $sequence_hash->{$prot_seq_id}{"taxonID"} = $leaf->taxon_id;
        # Taxon name
        $sequence_hash->{$prot_seq_id}{"taxon_name"} = $species_name;
        $species_name =~ s/\s/_/;
        #my $species_image_file = "http://localhost:3000/static/images/species_pictures/species_files/thumb_".$species_name.".png";
        my $species_image_file = "thumb_".$species_name.".png";
        # Taxon image
        $sequence_hash->{$prot_seq_id}{"taxon_image"} = $species_image_file;
	    # Gene description
        $sequence_hash->{$prot_seq_id}{"description"} = $gene4prot->description ;

        #push(@tmp_array, "<img src=''$species_image_file'' /> ".($leaf->taxon)->binomial);
        #push(@tmp_array,$leaf->stable_id );
        #push(@tmp_array,(defined $leaf->display_label ? $leaf->display_label: 'NaN'));
		#if(defined $leaf->display_label){ 
		#$symbol4family = $leaf->display_label;
		#	# grep description
		#	$symbol4family =~ s/-00\d//;
		#	my $grep_cmd = "grep -w \"$symbol4family\" /Users/fs9/Downloads/mart_export.txt";
		#	#print "$grep_cmd\n";
		#	my $grepped_line;# = `$grep_cmd`;
		#	chomp($grepped_line);
		#	my ($EnsemblGeneID,$EnsemblProteinID,$HGNCsymbol,$PFAMID,$WikiGeneDescription,$GOTermAccession,$GOTermName,$GOTermDefinition,$GOdomain,$GOTermEvidenceCode,$RfamID) = split("\t",$grepped_line);
		#	#print "grepped_line: $grepped_line";
		#	#$c->stash->{treefam}{pfamID} = (defined $PFAMID ? $PFAMID : "No Pfam ID");
		#	#$c->stash->{treefam}{Wikigenedescription} = (defined $WikiGeneDescription ? $WikiGeneDescription : "No WikiGene description");
		#	#$c->stash->{treefam}{GOTermAccession} = (defined $GOTermAccession ? $GOTermAccession : "No GOTermAccession");
		#	#
		#}
		#push(@tmp_array,$leaf->description );
        #$c->log->debug( "Found line:  ".join(",",@tmp_array)."  members" ) if $c->debug;
        #push(@array_of_arrays,\@tmp_array);
    }
    #$c->{'sequence_array_json'} =  encode_json \@array_of_arrays;
   
    #$c->log->debug('Family::Tree::get_summary decoded all sequences') if $c->debug; 
    # Count species
    #$c->{'numSpecies'} = keys(%species_count);
    foreach(@$all_leaves){ 
        $species_count->{$_->taxon->name} += 1;
    }
    $c->{'numSpecies'} = keys(%{$species_count});
    $c->{'speciesHash'} = \%{$species_count};
    $c->{'numSequences'} = scalar(@{$all_leaves});
    $c->{'tree_num_leaves'} = scalar(@{$all_leaves});
    
    #$c->log->debug( "Found " . keys(%species_count) . " species" ) if $c->debug;
    return 1;
}


sub get_sequence_annotations{
	my ($arg_ref) = @_;
	my $sequences_hash = $arg_ref->{sequences_hash};
	my $json_entry = $arg_ref->{sequence_array_json};
    my $hgnc_file= $arg_ref->{hgnc_file};
    my $ext_counts = $arg_ref->{ext_counts};
    my $wikigene_file= $arg_ref->{wikigene_file};
    my $pfam_file = $arg_ref->{pfam_file};
    my $hmmer_scores_file= $arg_ref->{hmmer_scores_file};
    my $uniprot_file= $arg_ref->{uniprot_file};
    my $db_adaptor= $arg_ref->{db_adaptor};
    
    my $counter = 0;
    my @all_sequences_array;
    print "Getting annotations for ".keys(%{$sequences_hash})." sequences\n";
    foreach my $ensembl_prot_id(keys(%{$sequences_hash})){
        #print "looking at $ensembl_prot_id\n";
        #######################################
        # HGNC
        #######################################
        my $gene4prot = $sequences_hash->{$ensembl_prot_id}{"geneID"};
        $sequences_hash->{$ensembl_prot_id}{"hgnc_hits"} = &get_hgnc_hits({"id" => $gene4prot, "file_to_search" => $hgnc_file, "db_adaptor" => $db_adaptor});
        $sequences_hash->{$ensembl_prot_id}{"wikigene_hits"} = &get_wikigene_hits({"id" => $gene4prot, "file_to_search" => $wikigene_file, "db_adaptor" => $db_adaptor });
        #print "WIKIGENE all: ".$sequences_hash->{$ensembl_prot_id}{"wikigene_hits"}."\n"; 
        #my $description_hits = &get_description_hits($leaf);
        #my $go_hits = &get_go_hits($leaf);
        $sequences_hash->{$ensembl_prot_id}{"pfam_hits"} = &get_pfam_hits({"id" => $ensembl_prot_id, "file_to_search" => $pfam_file, "pfam_counts" => \%{$ext_counts->{"pfam_counts"}}, "db_adaptor" => $db_adaptor});
        my $member_id4prot = $sequences_hash->{$ensembl_prot_id}{"member_id"};
        $sequences_hash->{$ensembl_prot_id}{"hmmer_hits"} = &get_hmmer_hits({"id" => $member_id4prot, "file_to_search" => $hmmer_scores_file, "db_adaptor" => $db_adaptor});
        $sequences_hash->{$ensembl_prot_id}{"uniprot_hits"} = &get_uniprot_hits({"id" => $ensembl_prot_id, "file_to_search" => $uniprot_file, "db_adaptor" => $db_adaptor});
    } 
    foreach my $ensembl_prot_id(keys(%{$sequences_hash})){
        my @pfam_hits;
        my @seq_array;
        #print Dumper $sequences_hash;
        #exit;
        # Seq order is: taxon, prot_id,
        push(@seq_array, (exists $sequences_hash->{$ensembl_prot_id}{"taxon_name"} && $sequences_hash->{$ensembl_prot_id}{"taxon_name"} ne '')?$sequences_hash->{$ensembl_prot_id}{"taxon_name"}: "NaN");
        push(@seq_array, $ensembl_prot_id);
        push(@seq_array, (exists $sequences_hash->{$ensembl_prot_id}{"hmmer_hits"} && $sequences_hash->{$ensembl_prot_id}{"hmmer_hits"} ne '')?$sequences_hash->{$ensembl_prot_id}{"hmmer_hits"}: "No hits" );
        if(exists $sequences_hash->{$ensembl_prot_id}{"pfam_hits"}){
                print "$ensembl_prot_id has ".keys(%{$sequences_hash->{$ensembl_prot_id}{"pfam_hits"}})." Pfam hits. \n";
                foreach my $acc(keys(%{$sequences_hash->{$ensembl_prot_id}{"pfam_hits"}})){
                    print "$acc hit has ".$ext_counts->{'pfam_counts'}{$acc}." other sequences (total: ".keys(%{$sequences_hash}).")....";
                    my $ratio = $ext_counts->{"pfam_counts"}{$acc} / keys(%{$sequences_hash});
                    if ( $ext_counts->{"pfam_counts"}{$acc} && $ratio < 0.05){
                        print "SKIP  this entry (ratio: $ratio) \n";
                        print "delete this entry ";
                        delete($sequences_hash->{$ensembl_prot_id}{"pfam_hits"}{$acc});
                    }
                    else{
                        print "\n";
                        push(@pfam_hits, $sequences_hash->{$ensembl_prot_id}{"pfam_hits"}{$acc}{"description"});
                    }
                }
        }
        
        if(!scalar(@pfam_hits)){
            push(@seq_array, "No hits");
        }
        else{
            push(@seq_array, join(",",@pfam_hits));
        }
        #my $hgnc_hits_string;
        #if(keys(%{$sequences_hash->{$ensembl_prot_id}->{"hgnc_hits"}})){ 
        #    $hgnc_hits_string =  join(",",keys(%{$sequences_hash->{$ensembl_prot_id}{"hgnc_hits"}}));
        #}
        #else{
        #    $hgnc_hits_string = "NaN";
        #}
        #push(@seq_array, $hgnc_hits_string);
        if(!exists($sequences_hash->{$ensembl_prot_id}{"description"}) ||  !defined($sequences_hash->{$ensembl_prot_id}{"description"}) || $sequences_hash->{$ensembl_prot_id}{"description"} eq ''){
            print "no description for $ensembl_prot_id\n";
            push(@seq_array, "NaN");
        }
        else{
            push(@seq_array, $sequences_hash->{$ensembl_prot_id}{"description"});
        }
        push( @all_sequences_array,\@seq_array);
       #last;
    }
       #die "we stop here, ok?"; 
    #$c->{"hgnc"}  = join(",",@hgnc_hits_array);
    #$c->{"pfam"}  = join(",",keys(%pfam_hits));
    $$json_entry = encode_json \@all_sequences_array;

    return 1;
    }

sub get_family_annotations{
	my ($arg_ref) = @_;
	my $c = $arg_ref->{superhash};
	my $tree = $arg_ref->{treeObject};
	my $species_count = $arg_ref->{species_count};
	my $sequences_hash = $arg_ref->{sequences_hash};
	my $previous_treefam_file = $arg_ref->{previous_treefam_file};
  ## get previous treefam annotations;
    my ($acc,$symbol,$desc,$seed,$full,$date1,$date2,$type,$full_desc,$comment) = &get_previous_treefam_annotation({"id" => $tree->stable_id, "file_to_search" => $previous_treefam_file});
    $c->{fam_symbol} = $symbol;
    $c->{fam_description} = $desc;
    $c->{fam_n_seed} = $seed;
    $c->{fam_n_full} = $full;
    #$c->{fam_full_description} = $full_desc;
    print "Symbol: $symbol, desc: $desc, n_seed: $seed, n_full: $full\n";
    #return 1;
    ### will do some summarizing things here as well
    # get hgnc symbols
    my $hgnc_string;
    my %pfam_hits;
    my %hgnc_hits;
    my $total_number_of_seqs =  keys(%{$sequences_hash});
    foreach my $seq(%{$sequences_hash}){
       # wikigene
       my %wikigene_hits;
        use Data::Dumper;
        #print Dumper $sequences_hash->{$seq}{'wikigene_hits'};
        print "WIKIGENE: ".keys(%{$sequences_hash->{$seq}{'wikigene_hits'}})." entries\n";
        if(keys(%{$sequences_hash->{$seq}{'wikigene_hits'}})){
            foreach my $wikigene_hit(keys(%{$sequences_hash->{$seq}->{wikigene_hits}})){
                $wikigene_hits{$wikigene_hit} = 1;
            print "WIKIGENE: add $wikigene_hit\n";
            }
        }
        # HGNC
        if(keys(%{$sequences_hash->{$seq}{'hgnc_hits'}})){
            foreach my $hgnc_hit(keys(%{$sequences_hash->{$seq}{hgnc_hits}})){
                $hgnc_hits{$hgnc_hit} = $sequences_hash->{$seq}{hgnc_hits}{$hgnc_hit};
            }
        }
        # PFAM
        if(keys(%{$sequences_hash->{$seq}{'pfam_hits'}})){
            foreach my $pfam_hit(keys(%{$sequences_hash->{$seq}{pfam_hits}})){
                $pfam_hits{$pfam_hit}{count}++;
                $pfam_hits{$pfam_hit}{id} = $sequences_hash->{$seq}{pfam_hits}{$pfam_hit}{hmm_name};
                $pfam_hits{$pfam_hit}{name} = $sequences_hash->{$seq}{pfam_hits}{$pfam_hit}{description};
            }
        }
        my @temp_hgnc_array;
        foreach my $id(keys(%hgnc_hits)){
            if(!exists $hgnc_hits{$id}){
                print "PARSE_ERROR: HGNC no id for $id\n";
            }
            else{
                push(@temp_hgnc_array, $id."->NaN");
            }
        }
        $c->{'hgnc'} = join(" ", @temp_hgnc_array);
        # PFAM
        my @temp_pfam_array;
        my @sorted_pfam_domains = sort { $pfam_hits{$a} <=> $pfam_hits{$b} } keys %pfam_hits;
        foreach my $id(@sorted_pfam_domains){
            my $percent = (int(($pfam_hits{$id}{count} * 100)/$total_number_of_seqs)) || 1;   
            push(@temp_pfam_array, $pfam_hits{$id}{name}."->".$pfam_hits{$id}{id}."->".$percent);
        }

        $c->{'pfam'} = join(" ",@temp_pfam_array);
        $c->{'wikigene'} = join(" ",keys(%wikigene_hits));
    }
   
    ## get overrepresented taxa
    my $total_no_seq_species = 0;

    my @species_array;
    my $top = 3;
    my $count_top = 1;
    foreach my $key (sort { $species_count->{$b} <=> $species_count->{$a}} keys %$species_count ){
           push(@species_array, $key."->".$species_count->{$key});
           print "TAXA_COUNT: ".$key."->".$species_count->{$key}."\n";
           last if $count_top++ >= $top;
    }
    $c->{'taxa_count'} = join(" ",@species_array);
    
    print "TAXA_COUNT: ".$c->{'taxa_count'};
    print "HGNC: ".$c->{'hgnc'};
    print "PFAM ".$c->{'pfam'};
   #exit; 
    return 1;
}


sub get_tree_nhx{
	my ($arg_ref) = @_;
	my $c = $arg_ref->{superhash};
	my $tree = $arg_ref->{treeObject};
	my $entry = $arg_ref->{entry};
    my $root_of_tree = $tree->root;
    $c->{'treeNhx'} = $root_of_tree->nhx_format;
  
    return defined($c->{'treeNhx'})? 1 : undef ;
}

sub get_alignment{
	my ($arg_ref) = @_;
	my $c = $arg_ref->{superhash};
	my $tree = $arg_ref->{treeObject};
	my $entry = $arg_ref->{entry};

# Get Alignment
    my $simpleAlign = $tree->get_SimpleAlign();
    my $fasta_string;
    foreach my $seq ( $simpleAlign->each_seq ){
                $fasta_string .= ">" . $seq->id . "\n" . $seq->seq . "\n";
    }
    $c->{'fasta_aln'} = $fasta_string;
    
  return 1;
  
}
sub get_hgnc_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
	my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select hgnc_id,app_symbol from hgnc_mapping where ensembl_gene_id = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my @hits;
 
    my %hgnc_hits;
### Parse
	while ( my ($hgnc_id,$app_symbol) = $extID2seq_sth->fetchrow_array() ){
        $hgnc_hits{$app_symbol} = $hgnc_id;
    }
    if(keys(%hgnc_hits)){
         print Dumper %hgnc_hits;
        #die "finished for now";
    } 
    return  \%hgnc_hits;
}

sub get_wikigene_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
    my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select wikigene_name,wikigene_id,wikigene_description from wikigene where ensembl_gene_id = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
 
    my %wikigene_hits;
    ### Parse
	while ( my ($wikigene_name,$wikigene_id,$wikigene_description) = $extID2seq_sth->fetchrow_array() ){
            $wikigene_hits{$wikigene_name}{id} = $wikigene_id;
            $wikigene_hits{$wikigene_name}{description} = $wikigene_description;
    } 
    return  \%wikigene_hits;
}

sub get_pfam_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
	my $pfam_counts = $arg_ref->{pfam_counts};
    my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select * from pfam_hits where seq_id = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
 
    my %wikigene_hits;
    ### Parse
    my %pfam_hits;
	while ( my ($ensID,$startA,$endA,$startB,$endB,$acc,$name,$type,$o1,$o2,$o3,$bitscore,$evalue,$o4,$clan) = $extID2seq_sth->fetchrow_array() ){
    #my $grep_command = "grep \"$id\" $file_to_search";
    #my @result_lines = `$grep_command`;
    #my @hgnc_hits;
### Parse
    #foreach (@result_lines){
        #print "fetched $_\n";    
        #my @array = split(/\t/,$_);
        #my ($ensID,$startA,$endA,$startB,$endB,$acc,$name,$type,$o1,$o2,$o3,$bitscore,$evalue,$o4,$clan) = split(/\s+/,$_);
        #print "$ensID,$startA,$endA,$startB,$endB,$acc,$name,$type,$o1,$o2,$o3,$bitscore,$evalue,$o4,$clan\n";
        $pfam_hits{$acc}{alignment_start} = $startA;
        $pfam_hits{$acc}{alignment_end} = $endA;
        $pfam_hits{$acc}{count}++;
        $pfam_hits{$acc}{description} = $name;
        $pfam_hits{$acc}{hmm_name} = $acc;
        $pfam_hits{$acc}{bitscore} = $bitscore;
        $pfam_hits{$acc}{confidence} = $evalue;
        $pfam_hits{$acc}{evalue} = $evalue;
        $pfam_counts->{$acc}++;
    } 
    #print Dumper $pfam_hits_href;
    #exit;
    return \%pfam_hits;
}

sub get_uniprot_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
    my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select * from uniprot_mapping where seq_id = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my %uniprot_hits;
### Parse
	while ( my ($seq_id,$db,$ext_seq_id) = $extID2seq_sth->fetchrow_array() ){
        $uniprot_hits{$db} = $ext_seq_id;
    } 
    return \%uniprot_hits;
}

sub get_hmmer_hits {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
    my $db_adaptor = $arg_ref->{db_adaptor};
    my $extID2seq_sth_ext = $db_adaptor->prepare('select evalue from hmmer_scores where member_id = ?');
	$extID2seq_sth_ext->bind_param(1,$id);
	my $extID2seq_sth = $extID2seq_sth_ext ;
    my $hmm_score; 
	$extID2seq_sth->execute() or die "SQL Error: $DBI::errstr\n";
    my %uniprot_hits;
### Parse
	while ( my ($evalue) = $extID2seq_sth->fetchrow_array() ){
        $hmm_score = $evalue;
        last;
    } 
    return $hmm_score;
}

sub get_previous_treefam_annotation {
	my ($arg_ref) = @_;
	my $id = $arg_ref->{id};
	my $file_to_search = $arg_ref->{file_to_search};
    my $grep_command = "grep \"^$id\" $file_to_search";
    my @result_lines = `$grep_command`;
    print "grep $grep_command\n";
### Parse
    my ($acc,$symbol,$desc,$empty,$seed,$full,$empty2,$date1,$date2,$type,$full_desc,$comment) = split(/\t/,$result_lines[0]);
    #print "$acc,$symbol,$desc,$empty,$seed,$full,$empty2,$date1,$date2,$type,full:$full_desc,$comment\n";
    #exit;
    return ($acc,$symbol,$desc,$seed,$full,$date1,$date2,$type,$full_desc,$comment);
}

sub read_species_info_table{
	my ($arg_ref) = @_;
	my $species_file = $arg_ref->{species_file};
	my $taxonIDMappings_href = $arg_ref->{taxonIDMappings_href};
	my @lines_of_file = `cat $species_file`;
	foreach my $species(@lines_of_file){
		chomp($species);
		my ($taxID,$sciName,$comName) = split(/\t/,$species);
		$taxonIDMappings_href->{$taxID}{'scientific_name'} = $sciName;
		$taxonIDMappings_href->{$taxID}{'common_name'} = $comName;
	}
	print "Found ".scalar(@lines_of_file)." species infos\n";
}
sub read_domain_info_file{
	my ($domain_file,$seqIDMappings_href) = (@_);
	my @lines_of_file = `cat $domain_file`;
	my %seq_counter;
	foreach my $line(@lines_of_file){
		chomp($line);
		my ($seqID,$alignment_start,$alignment_end,$envelope_start,$envelope_end,$hmm_acc,$hmm_name,$type,$hmm_start,$hmm_end,$hmm_length,$bit_score,$E_value,$significance,$clan,$predicted_active_site_residues) = split(/\s+/,$line);
		my $counter = (exists $seq_counter{$seqID})? $seq_counter{$seqID}+1 : 1;
		$seqIDMappings_href->{$seqID}{$counter}{"alignment_start"} = $alignment_start;
		$seqIDMappings_href->{$seqID}{$counter}{"alignment_end"} = $alignment_end;
		$seqIDMappings_href->{$seqID}{$counter}{"hmm_acc"} = $hmm_acc;
		$seqIDMappings_href->{$seqID}{$counter}{"hmm_name"} = $hmm_name;
		$seqIDMappings_href->{$seqID}{$counter}{"confidence"} = $E_value;
		$seq_counter{$seqID}++;
	}
	#print "Found ".scalar(@lines_of_file)." species infos\n";
	#print Dumper $seqIDMappings_href;
#exit;
}

sub read_sequence_length{
	my ($arg_ref) = @_;
	my $registry = $arg_ref->{registry};
	my $id_aref = $arg_ref->{terminal_ids};
	my $seqIDLengthMappings_href = $arg_ref->{seqIDLength};
    
    my $member_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'Member' );
	foreach my $id(@{$id_aref}){
    	print "READ_SEQUENCE_LENGTH: ENSEMBLPEP -> $id\n";
	    my $member = $member_adaptor->fetch_by_source_stable_id("ENSEMBLPEP",$id);
		my $member_length = $member->seq_length;
		$seqIDLengthMappings_href->{$id} = $member_length;
	}
	print "Found seq length for ".keys(%{$seqIDLengthMappings_href})." ids\n";
}

sub get_tree_phyloxml{
	my ($arg_ref) = @_;
	my $c = $arg_ref->{superhash};
	my $sequences_hash = $arg_ref->{sequences_hash};
	my $tree = $arg_ref->{treeObject};
	my $newickTree = $arg_ref->{treeNhx};
	my $registry = $arg_ref->{registry};
	my $entry = $arg_ref->{entry};
	my $treefam_id = $arg_ref->{treefam_name};
    my $phyloxml_file = $entry."test.phyloxml";
    #my $image_prefix = "http://static.ensembl.org/i/species/48/";
    my $image_prefix = "/static/images/species_pictures/species_files/";
    my %internalNodesMapping;
    my %leafNodesMapping;
    my %alreadyMapped;
    my $debug = 1;
    my ($fh, $tree_file) = tempfile( );
   
    # save newick string as file
    open my $nw_tree_out, ">", $tree_file or die "Could not open $tree_file\n";
    print {$nw_tree_out} $newickTree."\n";
    close $nw_tree_out || die "Could not close $tree_file\n";
    if(!-e $tree_file || ! -s $tree_file){
        warn "Problem saving tree to file $tree_file\n";
        return 0;
    } 

my $temp_tree = $treefam_id."temp_tree.nhx";
my $in = Bio::TreeIO->new( '-format' => 'nhx', '-file' => $tree_file );
while( my $tree = $in->next_tree ) {
	my @nodes = $tree->get_nodes();
	#print "Found ".scalar(@nodes)." nodes in total\n";
	die "Could not get nodes from tree" if !scalar(@nodes);
	foreach my $node(@nodes){
		
		if($node->is_Leaf){
			my $name = $node->id;
			my @tags = $node->get_all_tags;
			foreach(@tags){
				#print "get value for $_\n";
				my @values = $node->get_tag_values($_);
				$leafNodesMapping{$name}{$_}= $values[0];	
			}		
			#print "Found leaf: $name with \n";
			#print Dumper @tags;
			#exit;
		}
		# internal node
		# need to replace ids
		else{
			my $name = $node->id;
			next if !defined $name || $name eq "";
			my @tags = $node->get_all_tags;
			my $newID = $name."_".$alreadyMapped{$name}++;
			foreach(@tags){
				#print "get value for $_\n";
				my @values = $node->get_tag_values($_);
				$internalNodesMapping{$newID}{$_}= $values[0];	
			}	
            print "BOOTSTRAP saved: ".$node->bootstrap."\n";
			$internalNodesMapping{$newID}{'bootstrap'} = $node->bootstrap;
			$node->id($newID);
				
	
		}
		
	}
   #my $bio_phylo_tree = Bio::Phylo::Forest::Tree->new_from_bioperl($tr);

 #print Dumper %leafNodesMapping;
 #print Dumper %internalNodesMapping;
	my $out = new Bio::TreeIO(-file => ">$temp_tree", -format => 'nhx');
    $out->write_tree($tree);
#print $tree->to_string;
}
################################################################################
########     Parse as Bio::Phylo::Project
################################################################################
# we parse the newick as a project so that we end
# up with an object that has both the tree and the
# annotated OTUs
my $proj = parse(
    '-format' => 'newick',
    '-file' => $temp_tree,
    '-as_project' => 1,
);
#unlink($temp_tree);
# here we make the OTUs
my ($forest) = @{ $proj->get_items(_FOREST_) };
#print "Found ".scalar(@{ $proj->get_items(_FOREST_) })." nodes\n";
#my @nodes = @{ $proj->get_items(_NODE_) };
# it's easier to make a factory object for creating the annotations
my $fac = Bio::Phylo::Factory->new;
my $tree = $forest->first;

################################################################################
########     Species info
################################################################################

my $species_info_file = "species_info.table";
my %taxonIDMappings;
&read_species_info_table({"species_file" => $species_info_file, "taxonIDMappings_href" => \%taxonIDMappings});
die "problem reading species infos\n" if !keys(%taxonIDMappings);

################################################################################
########     Sequence length info
################################################################################
# Collect terminal taxa and get sequence length
my @all_terminal_taxa = @{$tree->get_terminals()};
my @all_terminal_ids = map {$_->get_name } @all_terminal_taxa;
my %seqIDLength;

&read_sequence_length({"registry" => $registry, "terminal_ids" => \@all_terminal_ids, "seqIDLength" => \%seqIDLength});

die "problem reading seq_length infos\n" if !keys(%seqIDLength);
#print Dumper %seqIDLength;
$tree->visit_depth_first(
		    '-in' => sub {
				    my $current_node = shift;			
				    return if !defined($current_node->get_name) || $current_node->get_name eq "";
				    print "looking at node ".$current_node->get_name."\n" if $debug;
				    ## INNER NODES
				    if($current_node->is_internal){ 
							    my $node_name = $current_node->id;
							    print "\tis an internal node\n" if $debug;
							    if(!exists $internalNodesMapping{$node_name}){
									    warn "Could not find ".$node_name." \n";
							    }
							    if(exists $internalNodesMapping{$node_name}{'bootstrap'}){
									    #print "setting bootstrap for $node_name\n";
									    #print "before: ".$current_node->get_generic('bootstrap');
									    #$current_node->set_generic( 'bootstrap' => $internalNodesMapping{$node_name}{'bootstrap'} );
									    #$current_node->set_score( $internalNodesMapping{$node_name}{'bootstrap'} );
									    #print "after: ".$current_node->get_generic('bootstrap')."\n";
									    #$current_node->score($internalNodesMapping{$node_name}{'bootstrap'});
									    my $bootstrap =$internalNodesMapping{$node_name}{'bootstrap'} ;
                                        if(!$bootstrap || $bootstrap eq ""){
                                            $bootstrap = 0;
                                        }
                                        print "setting bootstrap: $bootstrap\n"; 
                                        #$current_node->add_meta(
									        #$fac->create_meta(
									    	    #'-namespaces' => { 'pxml' => 'http://www.phyloxml.org/1.10/terms#' },
									    	    #'-triple' => { 'pxml:confidence' => $bootstrap },
									        #)
        							    #);
							    }
							    #return;
						        ## label events	
							    if(exists $internalNodesMapping{$node_name}{'D'}){
								    my $type_of_event =$internalNodesMapping{$node_name}{'D'} ;
								    #$current_node->score($internalNodesMapping{$node_name}{'bootstrap'});
								    my $event = _create_dummy_event(($type_of_event eq 'N')?1:0);
        						    print "EVENT: $event\n";
                                    #$event = '<speci><>';
                                    $current_node->add_meta(
                					    $fac->create_meta(
                        				    '-namespaces' => { 'pxml' => 'http://www.phyloxml.org/1.10/terms#' },
                        				    '-triple' => { 'pxml:events' => $event },
                					    )
        						    );
							    }
							    # replace name
							    $node_name =~ s/_\d+//g;
							    $current_node->set_name($node_name);
				    }
				    ## TERMINAL
				    if($current_node->is_terminal){
					    my $node_name = $current_node->id;
                        # we could rename leaves to a uniprot-style: CCNB1_HUMAN
					    print "\tterminal node: ".$node_name."\n";
                        if(exists $sequences_hash->{$node_name}{'uniprot_hits'}->{'UniProtKB-ID'}){
					        #if(exists $ens2uniprot{$node_name}){
                            print "\twe can replace it with ".$sequences_hash->{$node_name}{'uniprot_hits'}->{'UniProtKB-ID'}."!\n";
                            $current_node->id($sequences_hash->{$node_name}{'uniprot_hits'}->{'UniProtKB-ID'});
					    }
				    #print "\tis an terminal node\n";
							    if(exists $leafNodesMapping{$node_name}{'T'}){
								        my $tax_id = $leafNodesMapping{$node_name}{'T'};
								        my $scientific_name = $taxonIDMappings{$tax_id}{'scientific_name'};
								        my $common_name = $taxonIDMappings{$tax_id}{'common_name'};
	 							        if($tax_id eq "" || $scientific_name eq "" || $common_name eq ""){
										        warn "\tCould not find tax_id: $tax_id and/or scientific_name: $scientific_name and/or common_name: $common_name\n";
								        }
								        $scientific_name = "NaN" if $scientific_name eq "";
								        $tax_id = "0" if $tax_id eq "";
								        $common_name = "NaN" if $common_name eq "";
								        print "taxonomy with $tax_id,$scientific_name,$common_name\n";
								        my $image_file = "thumb_";
								        my $abused_scientific_name = $scientific_name;
								        $abused_scientific_name =~ s/ /_/g; 
								        my $image_path = "http://web-treefamdev.internal.sanger.ac.uk:3000/".$image_prefix."/thumb_".$abused_scientific_name.".png";
								        print "image_path: $image_path\n";
								        my $tax = _create_dummy_taxonomy($tax_id,$scientific_name,$common_name, $image_path);
								        $current_node->add_meta(
                				        $fac->create_meta(
                        			        '-namespaces' => { 'pxml' => 'http://www.phyloxml.org/1.10/terms#' },
                        			        '-triple' => { 'pxml:taxonomy' => $tax },
                					        )
        						        );
							    }
						        print "checking if $node_name is in seq hash ....";
                                if(exists $sequences_hash->{$node_name}{'pfam_hits'}){
									        print "yes\n";
                                            my $number_of_predictions = keys(%{$sequences_hash->{$node_name}{'pfam_hits'}});
										        print "\tfound $number_of_predictions domain predictions for $node_name\n";
									        if($number_of_predictions){
										        if(!exists $seqIDLength{$node_name}){ die "no seq length for $node_name\n";}
										        my $sequence_length = $seqIDLength{$node_name};
										        if($sequence_length eq ""){ die "no given seq length for $node_name\n";}
                                                my $arch = _create_dummy_architecture($sequence_length,\%{$sequences_hash->{$node_name}{'pfam_hits'}});
											        $current_node->add_meta(
											        $fac->create_meta(
												        '-namespaces' => { 'pxml' => 'http://www.phyloxml.org/1.10/terms#' },
												        '-triple' => { 'pxml:sequence' => $arch },
											        )
										        );
									    }
                                        else{
                                           my $sequence_length = $seqIDLength{$node_name};
                                                print "Writing empty architecture for $node_name with seq length: $sequence_length\n";
                                                my $arch = _create_dummy_empty_architecture($sequence_length);
											        $current_node->add_meta(
											        $fac->create_meta(
												        '-namespaces' => { 'pxml' => 'http://www.phyloxml.org/1.10/terms#' },
												        '-triple' => { 'pxml:sequence' => $arch },
											        )
										        );
 
 
                                            
                                            }
						    }
                            else{
                                                                                print "no\n";
                                }
    
				    }
		    }
    );
    # now write the output
    #print unparse( '-format' => 'phyloxml', '-phylo' => $proj );
    
    open my $phyloxmlOUT_FH, ">", $phyloxml_file or die "Could not open $phyloxml_file";
    print {$phyloxmlOUT_FH} unparse( '-format' => 'phyloxml', '-phylo' => $proj );
    close $phyloxmlOUT_FH or die "Could not close $phyloxml_file\n";

    ### HACK!
    # unable to set type of confidence
    my $in_place_replacement = "perl -i -pe 's/<confidence/<confidence type=\"bootstrap\"/g' $phyloxml_file";
    print "REPLACE $in_place_replacement\n";
    `$in_place_replacement`;
    # perl  -pe 's/<confidence/<confidence type="bootstrap"/g' 911644test.phyloxml
    $c->{'treephyloxml'} =  `cat $phyloxml_file`;
    $c->{'treephyloxml'} =~ s/\n//g ;
    #unlink($phyloxml_file);
	return 1; 
}

sub _create_dummy_architecture {
	my ($sequence_length,$domain_info_href) = (@_);
	my $domain_string = "<domain_architecture length=\"$sequence_length\">";
	foreach my $domain_number(keys(%$domain_info_href)){
		my($from,$to,$confidence, $hmm_name) = ($domain_info_href->{$domain_number}{'alignment_start'},
												$domain_info_href->{$domain_number}{'alignment_end'},
												$domain_info_href->{$domain_number}{'confidence'}, 
												$domain_info_href->{$domain_number}{'hmm_name'} );
		$domain_string .= "<domain from=\"$from\" to=\"$to\" confidence=\"$confidence\">$hmm_name</domain>";
	}
		$domain_string .= "</domain_architecture>";
}
sub _create_dummy_empty_architecture {
	my ($sequence_length) = (@_);
	my $domain_string = "<domain_architecture length=\"$sequence_length\">";
		$domain_string .= "<domain from=\"1\" to=\"2\" confidence=\"1\">No domains</domain>";
	$domain_string .= "</domain_architecture>";
}
sub _create_dummy_event {
	my $speciation = shift;
	return ($speciation == 1 )? '<speciations>1</speciations>' :  '<duplications>1</duplications>' ;
}

sub _create_dummy_bootstrap {
	my $bootstrap = shift;
	return "<confidence type=\"bootstrap\">$bootstrap</confidence>";
}
sub _create_dummy_taxonomy {
	my ($tax_id,$scientific_name,$common_name,$image_path) = (@_);
	## make it nice
	if($common_name =~ m/kangaroo/){
		print "before $common_name\n";
		
	}
	$scientific_name =~ s/\'//;
	$common_name =~ s/\'//;
	print "after $common_name\n";
	return  
	"<scientific_name>$scientific_name</scientific_name>
	<common_name>$common_name</common_name>";
    #<uri>$image_path</uri>";
} 
1;
