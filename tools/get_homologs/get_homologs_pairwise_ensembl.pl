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


my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -user   => 'anonymous',
    -dbname => 'ensembl_compara_68',
    -host   => 'ensembldb.ensembl.org',
    -port   => '5306'
);

my $speciesA = "takifugu_rubripes";
my $speciesB = "rattus_norvegicus";

# Get genome_db_ids
my $genome_db_adaptor = $db->get_GenomeDBAdaptor;
my $gdb1 = $genome_db_adaptor->fetch_by_registry_name("$speciesA");
my $gdb2 = $genome_db_adaptor->fetch_by_registry_name("$speciesB");

if(!$gdb1){die "Could not find $speciesA in genome_db\n";}
if(!$gdb2){die "Could not find $speciesB in genome_db\n";}
##### fetch GT by member
my $start_time = time();
print "fetching something\n";
my $homology_adaptor = $db->get_HomologyAdaptor;
#my $all_homologies = $homology_adaptor->fetch_all_by_genome_pair(1,2);

# Get corresponding method link species set
my $all_homologies = $homology_adaptor->fetch_all_by_genome_pair($gdb1->dbID,$gdb2->dbID);
#my $all_homologies = $homology_adaptor->fetch_all_by_genome_pair(1,2);
#my $all_homologies = $homology_adaptor->fetch_all_by_MethodLinkSpeciesSet(21);
print "\tdone fetching\n";

print(time()-$start_time, " secs to grab ", scalar(@$all_homologies), " using API\n");
exit;
my $counter = 0;
## For each homology
foreach my $this_homology (@$all_homologies) {
  # print the description (type of homology) and the
        ## subtype (taxonomy level of the event: duplic. or speciation)
    my $start_time = time();
      #print $this_homology->description, " [", $this_homology->subtype, "]\n";
  
        ## print the members in this homology
          my $members = $this_homology->get_all_Members();
          my $orthology_line;
            foreach my $this_member (@$members) {
      		#print $this_member->source_name, " ", $this_member->stable_id, "\n";
			# (", $this_member->genome_db->name, ")\n";
            $orthology_line .= $this_member->stable_id."\t";
            }
            $orthology_line .= $this_homology->description."\t".$this_homology->subtype;
          print "$orthology_line\n";
        #print(time()-$start_time, " secs to grab ", scalar(@$all_homologies), " using API\n");
        print "Count: ".$counter++."\n";
        }


