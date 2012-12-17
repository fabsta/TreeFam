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
use strict;
 use JSON;
use warnings;
use Getopt::Long;

my $method_link_id = 21;
my $homology_text;
my $limit = 10;
my $outputdir = "/lustre/scratch109/sanger/fs9/treefam/homologs";
my $result = GetOptions ("id=s" => \$method_link_id);
if($method_link_id =~ /input/){
        $method_link_id =~ s/input\.//;
    }

use Bio::EnsEMBL::Registry;
my $registry = 'Bio::EnsEMBL::Registry';
Bio::EnsEMBL::Registry->load_all("../../production_treefam_reg_conf.pl");

### Get the gene adaptor for human
my $genomedb_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'GenomeDB' );
if(!$genomedb_adaptor){
	die "Could not load gene_adaptor\n";
}
my $homology_adaptor = $registry->get_adaptor( 'TreeFam', 'Compara', 'Homology' );
if(!$homology_adaptor){
	die "Could not load gene_adaptor\n";
}
my $genetreenode_adaptor= $registry->get_adaptor( 'TreeFam', 'Compara', 'GeneTreeNode' );
if(!$genetreenode_adaptor){
	die "Could not load gene_adaptor\n";
}

my $methodlink_adaptor= $registry->get_adaptor( 'TreeFam', 'Compara', 'MethodLinkSpeciesSet' );
if(!$methodlink_adaptor){
	die "Could not load gene_adaptor\n";
}
my $method_link_object = $methodlink_adaptor->fetch_by_dbID($method_link_id) ;
my $genome_dbs = $method_link_object->species_set_obj->genome_dbs();
my ($species1,$species2) = ($genome_dbs->[0],$genome_dbs->[1]);
my ($species1_name,$species2_name) = ($species1->taxon->name, $species2->taxon->name);
$species1_name =~ s/ /_/;
$species2_name =~ s/ /_/;
my $outputfile = "$outputdir/";
my @sorted_species = sort {$a cmp $b} ($species1_name,$species2_name);

#if($species1_name cmp $species2_name){
    #$outputfile .= "$species2_name-$species1_name.txt";
#}
#else{
    #$outputfile .= "$species1_name-$species2_name.txt";
#}
$outputfile .= "$sorted_species[0]-$sorted_species[1].txt";
print "saving to $outputfile\n";
#exit;
my $speciesA = "takifugu_rubripes";
my $speciesB = "rattus_norvegicus";

# Get genome_db_ids
#my $genome_db_adaptor = $db->get_GenomeDBAdaptor;
#my $gdb1 = $genome_db_adaptor->fetch_by_registry_name("$speciesA");
#my $gdb2 = $genome_db_adaptor->fetch_by_registry_name("$speciesB");

#if(!$gdb1){die "Could not find $speciesA in genome_db\n";}
#if(!$gdb2){die "Could not find $speciesB in genome_db\n";}
##### fetch GT by member
my $start_time = time();
print "fetching something\n";
#my $all_homologies = $homology_adaptor->fetch_all_by_genome_pair(1,2);

# Get corresponding method link species set
#my $all_homologies = $homology_adaptor->fetch_all_by_genome_pair($gdb1->dbID,$gdb2->dbID);
#my $all_homologies = $homology_adaptor->fetch_all_by_genome_pair(1,2);
my $all_homologies = $homology_adaptor->fetch_all_by_MethodLinkSpeciesSet($method_link_id);
print "\tdone fetching\n";

print(time()-$start_time, " secs to grab ", scalar(@$all_homologies), " using API\n");
#exit;
my $counter = 0;
## For each homology
foreach my $this_homology (@$all_homologies) {
  # print the description (type of homology) and the
        ## subtype (taxonomy level of the event: duplic. or speciation)
    my $start_time = time();
      #print $this_homology->description, " [", $this_homology->subtype, "]\n";
        my $tree_node_id = $this_homology->tree_node_id();
        if($tree_node_id){ 
                #print "has tree node id $tree_node_id\n";
                my $genetreeNode = $genetreenode_adaptor->fetch_node_by_node_id($tree_node_id);
                if($genetreeNode){
                        my $homologue_genes = $this_homology->get_all_Members();
                        my $hashref = $genetreeNode->get_tagvalue_hash;
                        #print Dumper $hashref;
                        my ($type,$taxon,$bootstrap) =  ($genetreeNode->get_value_for_tag('node_type'),$genetreeNode->get_value_for_tag('taxon_name'), $genetreeNode->get_value_for_tag('bootstrap'));
                        $bootstrap = (defined($bootstrap))? $bootstrap : "NaN";
                        $taxon = (defined($taxon))? $taxon : "NaN";
                        $type = (defined($type))? $type : "NaN";
                        my ($gene1, $gene2) = ($homologue_genes->[0],$homologue_genes->[1]);
                        #print $gene1->stable_id."\t".($gene1->taxon)->name."\t".$gene2->stable_id."\t".($gene2->taxon)->name."\t";
                        #print "taxon: $taxon type: $type bootstrap: $bootstrap\n";
                        #$homology_text .= $gene1->stable_id."\t".($gene1->taxon)->name."\t".$gene2->stable_id."\t".($gene2->taxon)->name."\t$taxon\t$type\t$bootstrap\n";
                        $homology_text .= $gene1->stable_id."\t".($gene1->taxon)->name."\t".$gene2->stable_id."\t".($gene2->taxon)->name."\t$taxon\t$type\n";
                    }

        }

        #last if $counter++ > $limit;
        next;
        }
        write_to_file({"file_name" => $outputfile,"text" => $homology_text});

sub write_to_file{
	#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $file_name = $arg_ref->{file_name};
	my $text      = $arg_ref->{text};
	### OPENING FILE
	open my $out, '>>', $file_name or die "Couldn't open '$file_name': $!";
	### Writing file
	print {$out} $text;
	### CLOSING FILE
	close $out or die "Couldn't close '$file_name': $!";
	return 1;
}

