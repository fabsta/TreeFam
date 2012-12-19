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

### Connect to DB
my $registry = 'Bio::EnsEMBL::Registry';
my $registry_file = '../../registry/production_treefam_reg_conf.pl';
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
        
        #---------------------------------------------------------------------------------------------------
# SUMMARY
#---------------------------------------------------------------------------------------------------

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
my $no_species =keys(%species_count) ;
print "Found $no_species species\n";

my $no_duplicated_species = 0;
my $total_number_of_species = 101;
# Iterate over species
foreach my $species (keys(%species_count)){
	if($species_count{$species} > 1){
		$no_duplicated_species++;
	}

}
print "Found ".$no_duplicated_species." duplicated species\n";
my $species_percentage = int($no_species * 100 / $total_number_of_species);
print "Found ".$species_percentage." percent of species\n";

if($species_percentage > 95 && $no_duplicated_species < 3){
	print "SCG: $entry\n";
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
    }
#die "looked at last one (maybe 1) family\n";
#}

# update corresponding gene_tree_root entry!

die "Finished analysis\n";


#inserting into gene_tree_root_tag


