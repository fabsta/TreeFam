#!/usr/bin/perl
use Getopt::Long;
use strict;
use warnings;
use JSON;
use TreeFam::Production_modules;
use File::Temp qw/ tempdir tempfile /;

my $entry;
my $counter;
my $switch = "all" ;
my $write_switch = "all" ;
my $write_db = 1;
my $debug = 0;
my $clean_old_entries  = 1;
my $result = GetOptions ("id=s" => \$entry,    # numeric
                      "switch=s"   => \$switch,      # string
                      "write_switch=s"  => \$write_switch,
                        "counter=s" => \$counter);

my $mappings_dir = "/nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/treefam_tools/make_wikipedia_mappings/";
my $names_and_entrez = "$mappings_dir/names_and_entrez.txt";
my $human_and_entrez = "$mappings_dir/results_entrez_homo_sapiens";

print "Check existence of mapping files" if $debug;
die "File desc missing\n" if(!-e $names_and_entrez|| !-s $names_and_entrez);
die "File go missing\n" if(!-e $human_and_entrez|| !-s $human_and_entrez);

print "Get DB connection\n" if $debug;

### Connect to DB
my $registry = 'Bio::EnsEMBL::Registry';
my $registry_file = '/nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/production_treefam_reg_conf.pl';
#$ENV{'ENSEMBL_REGISTRY'} = '/nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/production_treefam_reg_conf.pl';
TreeFam::Production_modules::get_db_connection({"registry" => $registry,"registry_file" => $registry_file});
die "Could not connect to DB. Check registry in $registry_file\n" if !$registry;
my $genetree_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'GeneTree' );
if(!$genetree_adaptor){
    die "Could not load genetree_adaptor\n";
}
my %superhash;


my $entrez2name_hash  = read_names_entrez($names_and_entrez);
die "problem reading $names_and_entrez file\n" if !keys(%{$entrez2name_hash});
my $gene2entrez_hash  = read_genes_entrez($human_and_entrez);
die "problem reading $human_and_entrez file\n" if !keys(%{$gene2entrez_hash});

print "read ".keys(%{$entrez2name_hash})." entrez entries and ".keys(%{$gene2entrez_hash})." gene2entrez entries\n";
exit;

print "ok, we have ".keys(%{$entrez2name_hash})." entrez2name and ".keys(%{$gene2entrez_hash})." gene2entrez entries\n" if $debug;
#my @families = ("11922233"); # TF101001
my @families = ("12619130","12255948", "11975910","12698788","12506762","11655023","12405244","11595034");
#TREE:
my @all_trees = @{$genetree_adaptor->fetch_all( -CLUSTERSET_ID => "default", -TREE_TYPE => "tree")};
#print "Found ".scalar(@all_trees)." trees\n";
my %names4family;
foreach my $gt (@all_trees){
    my $entry = $gt->root_id;
    my $treefam_name = $gt->stable_id;
    my %sequences_hash;	
    my %species_count;	
    my $mark_family = 0;
    die "wait, entry is invalid: $entry\n" if !defined($entry) || $entry eq '';
	print "Looking at: $entry\n" if $debug;
    print "\tGet Tree object\n" if $debug;
    # get all sequences
    my $all_leaves = ( $gt )->get_all_leaves();
    print "we have ".scalar(@{$all_leaves})." sequences to look at\n" if $debug;
    foreach my $sequence (@{$all_leaves}){
        my $gene4prot = $sequence->gene_member;
        my $gene_id = $gene4prot->stable_id;
        next if $gene4prot->taxon_id != 9606;
        print "\tsearching with gene_id $gene_id\n" if $debug;
        if(exists $gene2entrez_hash->{$gene_id}){
            my $entrez4gene =$gene2entrez_hash->{$gene_id} ;
            if(exists $entrez2name_hash->{$entrez4gene}){
                my $name4entrez =$entrez2name_hash->{$entrez4gene}  ;
                print "gene_id: $gene_id entrez: $entrez4gene name: $name4entrez\n" if $debug; 
                print "$treefam_name\t$name4entrez\n"; 
                #$names4family{$treefam_name}{$name4entrez} = 1;
                #last;
            }     
        }
    }
} 



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
sub read_genes_entrez{
	my ($file) = (@_);
	my %gene2entrez;
    my $counter = 0;
    my @lines_of_file = `cat $file`;
	foreach my $line(@lines_of_file){
		chomp($line);
		my ($id,$gene,$entrez,$wiki,$desc) = split(/\s+/,$line);
        $gene =~ s/\"//g;
	    $entrez =~ s/\"//g;
        next if $entrez eq "NA";

	    #print  "$id,$gene,$entrez,$wiki,$desc\n";
        #exit if $counter++ > 3;
        $gene2entrez{$gene} = $entrez;
    }
	print "Found ".keys(%gene2entrez)." entrez infos\n";
    return \%gene2entrez;
}
sub read_names_entrez{
	my ($file) = (@_);
	my %entrez2name;
    my @lines_of_file = `cat $file`;
	foreach my $line(@lines_of_file){
		chomp($line);
		my ($wikiname,$entrez) = split(/\t/,$line);
	    $entrez2name{$entrez} = $wikiname;
    }
	print "Found ".keys(%entrez2name)." entrez infos\n";
    return \%entrez2name;
}
