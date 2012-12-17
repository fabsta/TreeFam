#!/usr/bin/perl
use Data::Dumper;
use Scalar::Util;
use Getopt::Long;
use strict;
use warnings;
use JSON;
use File::Temp qw/ tempdir tempfile /;
use TreeFam::Production_modules;
use JSON;
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

#print "Get DB connection\n";


### Connect to DB
my $registry = 'Bio::EnsEMBL::Registry';
my $registry_file = '/nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/production_treefam_reg_conf.pl';
#$ENV{'ENSEMBL_REGISTRY'} = '/nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/production_treefam_reg_conf.pl';
TreeFam::Production_modules::get_db_connection({"registry" => $registry,"registry_file" => $registry_file});
die "Could not connect to DB. Check registry in $registry_file\n" if !$registry;


my $genome_db_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'GenomeDB' );
if(!$genome_db_adaptor){
    die "Could not load genome_db_adaptor\n";
}
my $member_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'Member' );
if(!$member_adaptor){
    die "Could not load member_adaptor\n";
}

my $all_genomes = $genome_db_adaptor->fetch_all();
my @all_genomes_array;
foreach my $genome (@{$all_genomes}){
    my $genome_hash;
    #print "looking at ".$genome->name."\n";
    my ($id,$name,$common_name,$classification ) =  ($genome->dbID,$genome->name, ($genome->taxon)->common_name, ($genome->taxon)->classification);
    if((!$name || $name eq "") || (!$common_name || $common_name eq "") || (!$classification || $classification eq "")){
        #print "There is a problem with $name,$common_name,$classification \n";
    }
    my @class = split(" ", $classification);
    my @rev_class = reverse(@class);
    
    my $new_class = join(';',@rev_class[0..5]);
    my ($first,$second) = split("_",$name);

    $genome_hash->{"name"} = ucfirst($first)."_".$second;
    $genome_hash->{"id"} = $id;
    #$genome_hash->{"common_name"} = $common_name;
    $genome_hash->{"common_name"} = (!$common_name || $common_name eq "")? "NaN" : $common_name ;
    $genome_hash->{"classification"} = $new_class;
    #my $no_of_sequences = $member_adaptor
    push(@all_genomes_array, $genome_hash);
}

#print Dumper @all_genomes_array;
my $json_string = encode_json(\@all_genomes_array);
print $json_string;

