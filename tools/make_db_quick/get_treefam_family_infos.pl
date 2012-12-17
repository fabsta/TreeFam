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
 use JSON;
use warnings;

#my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
#                                                 -user   => 'root',
#                                                 -dbname => 'treefam_homology_67x',
#                                                 -host   => 'localhost',
#                                                 -pass   => '123',
#                                                 -port   => '3306'
#                                               );

my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(-user   => 'treefam_admin',
-dbname => 'treefam_homology_67hmm',
-host   => 'web-mei-treefam',
-pass => 'treefam_king1982',
-port => '3365');
#my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
#    -user   => 'anonymous',
#    -dbname => 'ensembl_compara_67',
#    -host   => 'ensembldb.ensembl.org',
#    -port   => '5306'
#);

my $genetree_adaptor = $db->get_GeneTreeAdaptor;

my $tree = $genetree_adaptor->fetch_by_stable_id("TF101001");

if ( !defined($tree) )
{
    print "Could not get tree\n";
}
else
{
    print "get root\n";
    my $root_of_tree = $tree->root;
}

#my @tags = $tree->get_all_tags();
#print "got ".scalar(@tags)." tags\n";
#die;

#$ns_node->_load_tags();
#my $root_of_tree = $tree->root;

#my $tagvalue_hashref = $tree->get_tagvalue_hash();
#if ( !keys(%$tagvalue_hashref) )
#{
#    die "Could not get tagvalue_hashref for tree\n";
#}
my $summaryData;
my $c;


	#$summaryData->{numSequences}              = $tree->get_value_for_tag('gene_count');
	
   #die; 
    
    $summaryData->{tree_max_branch}           = $tree->get_value_for_tag('tree_max_branch');
    $summaryData->{tree_num_human_peps}       = $tree->get_value_for_tag('tree_num_human_peps');
    $summaryData->{tree_num_dup_nodes}        = $tree->get_value_for_tag('tree_num_dup_nodes');
    $summaryData->{aln_method}                = $tree->get_value_for_tag('aln_method');
    
    $summaryData->{aln_percent_identity}      = $tree->get_value_for_tag('aln_percent_identity');
    $summaryData->{aln_num_residues}          = $tree->get_value_for_tag('aln_num_residues');
    $summaryData->{tree_num_spec_nodes}       = $tree->get_value_for_tag('tree_num_spec_node');
    $summaryData->{aln_length}                = $tree->get_value_for_tag('aln_length');
    $summaryData->{aln_runtime}               = $tree->get_value_for_tag('aln_runtime');
    $summaryData->{tree_max_length}           = $tree->get_value_for_tag('tree_max_length');
    $summaryData->{buildhmm_runtime_msec}     = $tree->get_value_for_tag('buildhmm_runtime_msec');
    $summaryData->{njtree_phyml_runtime_msec} = $tree->get_value_for_tag('njtree_phyml_runtime_msec');
    $summaryData->{orthotree_runtime_msec}    = $tree->get_value_for_tag('orthotree_runtime_msec');
    $summaryData->{tree_num_leaves}           = $tree->get_value_for_tag('tree_num_leaves');
    
    $summaryData->{tree_info}{seed}              = $tree->get_value_for_tag('treenhx');
    $summaryData->{sequence_array_json}          = $tree->get_value_for_tag('sequence_array_json');
	
	$summaryData->{numSpecies}                = $tree->get_value_for_tag('numspecies');
	$summaryData->{fasta_aln}                = $tree->get_value_for_tag('fasta_aln');
	$summaryData->{"homologs_array_json"} = $tree->get_value_for_tag('homologs_array_json'); 
	
	
	print Dumper $summaryData->{fasta_aln};
	print Dumper $summaryData->{tree_info}{seed};
	print Dumper $summaryData->{sequence_array_json};
	print Dumper $summaryData->{"homologs_array_json"};
	
	print "got everything!!!!\n";
	die "Finished\n";

