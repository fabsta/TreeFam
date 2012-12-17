#
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
package TreeFam::HomologyHelper;

use strict;
use warnings;
 
sub get_member_object{
	my ($db,$source_id) = (@_);
	my $member_object = $db->get_MemberAdaptor;
	my $member = $member_object->fetch_by_source_stable_id(undef,$source_id);
	print "getting member object\n";
	return defined($member)?$member:undef;
}
sub get_genetree_object{
	my ($db,$genetree_id) = (@_);
	my $genetree_object = $db->get_GeneTreeAdaptor;
	my $genetree = $genetree_object->fetch_by_stable_id($genetree_id);
	return defined($genetree)?$genetree:undef;
}
sub get_all_species_for_family{
	my ($db,$genetree_id) = (@_);
	my %have_species;
	my @all_species;
	my $gt = &get_genetree_object($db,$genetree_id);
	if(!$gt){return undef;}
	my $all_members = ($gt->root)->get_all_leaves();
	foreach my $member(@{$all_members}){
	my $species_name = $member->taxon->name;
		push(@all_species, $species_name)  if !exists $have_species{$species_name};
		$have_species{$species_name} = 1;
	}

    return scalar(@all_species)?\@all_species:undef;
}
sub get_all_homologs_for_family{
	my ($db,$source_genetree_node, $type) = (@_);
	my $homology_adaptor = $db->get_HomologyAdaptor;
	my $homologies = $homology_adaptor->fetch_all_by_tree_node_id($source_genetree_node);
	my $encoded_homologies = encode_homologies($homologies,$type);
	return defined($encoded_homologies)?$encoded_homologies:undef
}
sub get_member_pair_species{
	my ($db,$source_member_object, $second_species) = (@_);
	my $homology_adaptor = $db->get_HomologyAdaptor;
	my $homologies = $homology_adaptor->fetch_all_by_Member_paired_species($source_member_object, $second_species);
    return defined($homologies)?$homologies:undef
}
sub get_homologs_for_gene{
	my ($db,$source_member_object, $type, $homology_type) = (@_);
	my $memberID = $source_member_object->dbID;
	# the api requires member ids (int) rather than member objects
	my $homology_adaptor = $db->get_HomologyAdaptor;
	print "have homology adaptor";
	my $member_adaptor = $db->get_MemberAdaptor;
	#my $member = $member_adaptor->fetch_by_source_stable_id("ENSEMBLGENE", "ENSG00000229314");
	my $homologies = $homology_adaptor->fetch_all_by_Member($source_member_object);
	#my $homologies = $homology_adaptor->fetch_all_by_Member($source_member_object);
	#return undef if !scalar(@$homologies);
	my $genetree = &get_genetree_for_member($db,$source_member_object);
	my $tf_family = $genetree->[0]->stable_id."\n";
	
	my $encoded_homologies = encode_homologies($homologies,$type, $homology_type,$tf_family);
	#print Dumper $encoded_homologies;
	return defined($encoded_homologies)?$encoded_homologies:undef
}
sub get_homology_relation_for_gene_pair{
	my ($db,$source_member_object, $source_member2_object) = (@_);
	my ($memberID1,$memberID2) = ($source_member_object->dbID, $source_member2_object->dbID);
	# the api requires member ids (int) rather than member objects
	my $homology_adaptor = $db->get_HomologyAdaptor;
	my $homologies = $homology_adaptor->fetch_by_Member_id_Member_id($memberID1,$memberID2);
    return defined($homologies)?$homologies:undef;
}
sub get_genetree_for_member{
	my ($db,$member) = (@_);
	print "\tsearching for genetree \n";
	my $genetree_adaptor = $db->get_GeneTreeAdaptor;
	my $genetrees = $genetree_adaptor->fetch_all_by_Member($member, -CLUSTERSET_ID => "default");
	return scalar(@{$genetrees})? $genetrees : undef;
}	
sub get_homologs_for_family_level{
	my ($db,$gt_object, $taxon_level) = (@_);
	print "\tsearching for $taxon_level with ".$gt_object->stable_id."\n";
	my $root = $gt_object->root;
	## get all nodes that match taxon level 
	# this can be multiple
	#my $homologies = $source_genetree_object->get_all_nodes_by_tagvalue('taxon_name'=>$taxon_level);
	my $taxon_level_nodes= $root->get_all_nodes_by_tag_value('taxon_name'=>$taxon_level);
	my $homology_adaptor = $db->get_HomologyAdaptor;
	my @all_homologies;
	foreach my $node(@{$taxon_level_nodes}){
		my $homologies4node = $homology_adaptor->fetch_all_by_tree_node_id($root->node_id);
		push(@all_homologies, @{$homologies4node});
	}
    return defined(\@all_homologies)?\@all_homologies:undef
}
sub encode_homologies{
	my ($homologies,$type, $homology_type, $tf_family) = (@_);
	my @results_array;
	my @sequences_array;
	my %sequences_remember_hash;
	foreach my $this_homology (@$homologies) {
		if(defined($homology_type)  && $homology_type ne "" && $this_homology->description ne $homology_type){
			#print "skipping type: ".$this_homology->description."\n";
			next;
		}
		my $homologue_genes = $this_homology->gene_list;
    	if(scalar(@{$homologue_genes}) > 2){
			die "too many genes: ".join(" and ", @$homologue_genes), " are ", $this_homology->description, "\n";
		}
		else{
			my ($source,$target) = ($homologue_genes->[0],$homologue_genes->[1]);
			my %pair_hash;
			my ($source_protein_id,$target_protein_id) = ($source->get_canonical_Member->stable_id,$target->get_canonical_Member->stable_id);
			my ($source_protein_sequence,$target_protein_sequence) = ($source->get_canonical_Member->sequence(),$target->get_canonical_Member->sequence());
			# source
					$pair_hash{"source"} = {
							"protein_id" => $source_protein_id,
							"species" => $source->taxon->name,
							"id" => $source->stable_id,
							"sequence" => $source_protein_sequence
						};
			# target
					$pair_hash{"target"} = {
							"protein_id" => $target_protein_id,
							"species" => $target->taxon->name,
							"id" => $target->stable_id,
							"sequence" => $target_protein_sequence
						};
			# type
				$pair_hash{"tf_family"} = $tf_family;
				$pair_hash{"type"} = $this_homology->description();
				$pair_hash{"subtype"} = $this_homology->subtype();
				
				#my $target_fasta = ">".$target_protein_id."\n".$target_protein_sequence;	
				#my $source_fasta = ">".$source_protein_id."\n".$source_protein_sequence;
				#push(@sequences_array, \$target_fasta) if !exists $sequences_remember_hash{$target_protein_id};
				#push(@sequences_array, \$source_fasta) if !exists $sequences_remember_hash{$source_protein_id};
				push(@results_array, \%pair_hash);		
			# Remember we saved the sequences
			#$sequences_remember_hash{$target_protein_id} = 1;
			#$sequences_remember_hash{$source_protein_id} = 1;
		}
  	}
	#print "found ".scalar(@results_array)." homologies and ".scalar(@sequences_array)." sequences\n";
  #return \@results_array,\@sequences_array;
  return \@results_array;
}


sub get_homologs_for_family_orthoxml{
	my ($db,$source_genetree, $write_type) = (@_);
	my ($w,$string_handle,$file_handle);	
	my $file_name = "output.xml";
	if($write_type eq "file"){
		$file_handle = IO::File->new($file_name, 'w');
  		$w = Bio::EnsEMBL::Compara::Graph::OrthoXMLWriter->new(	-SOURCE => 'Ensembl', -SOURCE_VERSION => 67, -HANDLE => $file_handle);
	}
	elsif($write_type eq "string"){
		$string_handle = IO::String->new();
  		$w = Bio::EnsEMBL::Compara::Graph::OrthoXMLWriter->new(
    		-SOURCE => 'Ensembl', -SOURCE_VERSION => 67, -HANDLE => $string_handle
  		);
	}
	else{ die "Wrong type for writing given for get_homology_for_family_orthoxml: $write_type\n";}
  	$w->write_trees($source_genetree);
  	$w->finish(); #YOU MUST CALL THIS TO WRITE THE FINAL TAG
	if($write_type eq "string"){
  		my $xml_scalar_ref = $string_handle->string_ref();
		return defined($xml_scalar_ref)?$xml_scalar_ref:undef;	
	}
	elsif($write_type eq "file"){
		$file_handle->close();
		return (-e $file_name && -s $file_name)? $file_name:undef;
	}
	else{
		die "wrong option\n";
	}
}

sub check_existence {
	my ($object,$type) = shift;
	if(!$object){ die "Could not get member object for $type. Stopping\n"}
}
1;
