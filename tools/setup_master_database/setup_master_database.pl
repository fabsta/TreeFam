

## production run of the TreeFam pipeline
# set master database
#

#Testing
#Testing loading from fasta files
# database infos
use strict;
use warnings;
use DBI;
use TreeFam::Production_modules;



########################################################
my $write = 1;
my $write_orthologs = 1 ;
my $db_host = "web-mei-treefam";
my $db_user = $ENV{"TREEFAM_DB_USER"};
my $db_pass = $ENV{"TREEFAM_DB_PASS"};
my $db_port = "3365";
my $master_db_name = "treefam_master9";
die "No username provided\n" if $db_user eq '';
die "No password provided\n" if $db_pass eq '';
my $ENSEMBL_CVS_ROOT_DIR = $ENV{"ENSEMBL_CVS_ROOT_DIR"};
my $ensembl_version = '69';
my $ensembl_compara_db = "ensembl_compara_$ensembl_version";
my $ensembl_metazoa_compara_db = "ensembl_compara_metazoa_16_69";
my $ensembl_fungi_compara_db = "ensembl_compara_fungi_16_69";
my $ensembl_plants_compara_db = "ensembl_compara_plants_16_69";

my %fungi = (
"schizosaccharomyces_pombe" => 1,
);
my %plants = (
"arabidopsis_thaliana" => 1,
);
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
##########################################
my %all_species;
my $species_file = "species_file";
my $species_count = 0;
# Create master database
my $create_db_cmd = "mysql -h $db_host -P $db_port -u $db_user -p$db_pass -e 'create database $master_db_name'";
system($create_db_cmd) if $write;# or die "Could not create master db: $create_db_cmd\n";
print "$create_db_cmd\n";

# load database schema
my $load_db_schema_cmd = " mysql -h $db_host -P $db_port -u $db_user -p$db_pass -D $master_db_name < $ENSEMBL_CVS_ROOT_DIR/ensembl-compara/sql/table.sql";
system($load_db_schema_cmd) if $write ;# or die "Could not load schema for db: $load_db_schema_cmd\n";
print "$load_db_schema_cmd\n";

# Import Method link data
my $load_method_links_cmd = "mysqlimport --local -h $db_host -P $db_port -u $db_user -p$db_pass  $master_db_name $ENSEMBL_CVS_ROOT_DIR/ensembl-compara/sql/method_link.txt";
system($load_method_links_cmd) if $write;# or die "Could not load method link table: $load_method_links_cmd\n";
print "$load_method_links_cmd\n";
#Put all fasta files in directories

#Taxon names have to match NCBI taxonomy
#Functions to be changed

my $dump_ncbi_tables_cmd = "mysqldump -h ensembldb.ensembl.org --skip-lock-tables -P5306 -uanonymous $ensembl_compara_db ncbi_taxa_name ncbi_taxa_node | mysql -h $db_host -P $db_port -u $db_user -p$db_pass  $master_db_name ";
system($dump_ncbi_tables_cmd) if $write;# or die "Could not dump ncbi tables: $dump_ncbi_tables_cmd\n";
print "$dump_ncbi_tables_cmd\n";


########################################################
#Load ensembl species
########################################################
#!! make sure the production_reg_conf points to the correct database 
#     (changes in master db)
#                         or previous compara version (e.g. 67 -> 68)
#
#! also, make sure the offset is correctly set 
#Offset:
#     by default, the first newly added genome gets genome_id 1, next one 2, and so on. To be able to reuse information
#from a master database, we need to make sure genome_db_ids are the same 

#! manually, but we need to change this
# Check
#     New release of genomes?
my $update_genome_file = "$ENSEMBL_CVS_ROOT_DIR/ensembl-compara/scripts/pipeline/update_genome.pl";
my $registry_file = "$ENSEMBL_CVS_ROOT_DIR/treefam/scripts/pipeline/production_treefam_reg_conf.pl";

my $update_genome_script = "perl $update_genome_file  --reg_conf $registry_file --compara treefam_master --species ";
# add the species in each iteration
{
# First, read all genomes, we have already imported.i
#$all_genome_dbs = $genome_db_adaptor->fetch_all();
#my $dbh_ensembl = DBI->connect("dbi:mysql:$ensembl_compara_db:ensembldb.ensembl.org:5306", 'anonymous', '') or die "could not connect to ensembl database\n";
#my $query_genome_db = "select * from genome_db;";
#my $query_handle = $dbh_ensembl->prepare($query_genome_db);
 #$query_handle->execute();
#my ($genome_db_id, $taxon_id, $name, $assembly, $assembly_default, $genebuild, $locator);
#

# Now read all available species from the ensembl database;
my $dbh_ensembl = DBI->connect("dbi:mysql:$ensembl_compara_db:ensembldb.ensembl.org:5306", 'anonymous', '') or die "could not connect to ensembl database\n";
my $query_genome_db = "select * from genome_db;";
my $query_handle = $dbh_ensembl->prepare($query_genome_db);
 $query_handle->execute();
my ($genome_db_id, $taxon_id, $name, $assembly, $assembly_default, $genebuild, $locator);
# BIND TABLE COLUMNS TO VARIABLES
$query_handle->bind_columns(\$genome_db_id, \$taxon_id, \$name, \$assembly, \$assembly_default, \$genebuild, \$locator);

while($query_handle->fetch()) {
    next if $name eq 'ancestral_sequences';
    #print "$genome_db_id, $taxon_id, $name, $assembly, $assembly_default, $genebuild, $locator\n";
    next if exists $all_species{$name};
    print "$name\n";
    my $species_specific_cmd = $update_genome_script." \"$name\" ";
    print "$species_specific_cmd\n";
    system("$species_specific_cmd") if $write;
    $all_species{$name} = 1;
    $species_count++ if !exists $all_species{$name};;
#    $species_count_hash{'EnsEMBL'}{$name} = 1;
}
}

########################################################
#Load ensembl metazoa species
########################################################
# Now read all available species from the ensembl database;
{
my $dbh_ensembl_metazoa = DBI->connect("dbi:mysql:$ensembl_metazoa_compara_db:mysql.ebi.ac.uk:4157", 'anonymous', '') or die "could not connect to ensembl metazoa database\n";
my $query_genome_db = "select * from genome_db;";
my $query_handle = $dbh_ensembl_metazoa->prepare($query_genome_db);
 $query_handle->execute();
my ($genome_db_id, $taxon_id, $name, $assembly, $assembly_default, $genebuild, $locator);
# BIND TABLE COLUMNS TO VARIABLES
$query_handle->bind_columns(\$genome_db_id, \$taxon_id, \$name, \$assembly, \$assembly_default, \$genebuild, \$locator);

while($query_handle->fetch()) {
    next if $name eq 'ancestral_sequences';
    next if exists $all_species{$name};
    #print "$genome_db_id, $taxon_id, $name, $assembly, $assembly_default, $genebuild, $locator\n";
    print "$name\n";
    my $species_specific_cmd = $update_genome_script." \"$name\" ";
    print "$species_specific_cmd\n";
    system("$species_specific_cmd") if $write;
    $all_species{$name} = 1;
    $species_count++ if !exists $all_species{$name};;
#    $species_count_hash{'EnsEMBL'}{$name} = 1;
}
}
########################################################
#Load ensembl FUNGI species
########################################################
# Now read all available species from the ensembl database;
{
my $dbh_ensembl_fungi = DBI->connect("dbi:mysql:$ensembl_fungi_compara_db:mysql.ebi.ac.uk:4157", 'anonymous', '') or die "could not connect to ensembl fungi database\n";
my $query_genome_db = "select * from genome_db;";
my $query_handle = $dbh_ensembl_fungi->prepare($query_genome_db);
 $query_handle->execute();
my ($genome_db_id, $taxon_id, $name, $assembly, $assembly_default, $genebuild, $locator);
# BIND TABLE COLUMNS TO VARIABLES
$query_handle->bind_columns(\$genome_db_id, \$taxon_id, \$name, \$assembly, \$assembly_default, \$genebuild, \$locator);

while($query_handle->fetch()) {
    next if $name eq 'ancestral_sequences';
    next if exists $all_species{$name};
    next if !exists $fungi{$name};
    #print "$genome_db_id, $taxon_id, $name, $assembly, $assembly_default, $genebuild, $locator\n";
    print "$name\n";
    my $species_specific_cmd = $update_genome_script." \"$name\" ";
    print "$species_specific_cmd\n";
    system("$species_specific_cmd") if $write;
    $all_species{$name} = 1;
    $species_count++ if !exists $all_species{$name};;
#    $species_count_hash{'EnsEMBL'}{$name} = 1;
}
}

########################################################
#Load ensembl PLANTS species
########################################################
# Now read all available species from the ensembl database;
{
my $dbh_ensembl_plants = DBI->connect("dbi:mysql:$ensembl_plants_compara_db:mysql.ebi.ac.uk:4157", 'anonymous', '') or die "could not connect to ensembl plants database\n";
my $query_genome_db = "select * from genome_db;";
my $query_handle = $dbh_ensembl_plants->prepare($query_genome_db);
 $query_handle->execute();
my ($genome_db_id, $taxon_id, $name, $assembly, $assembly_default, $genebuild, $locator);
# BIND TABLE COLUMNS TO VARIABLES
$query_handle->bind_columns(\$genome_db_id, \$taxon_id, \$name, \$assembly, \$assembly_default, \$genebuild, \$locator);

while($query_handle->fetch()) {
    next if $name eq 'ancestral_sequences';
    next if exists $all_species{$name};
    next if !exists $plants{$name};
    #print "$genome_db_id, $taxon_id, $name, $assembly, $assembly_default, $genebuild, $locator\n";
    print "$name\n";
    my $species_specific_cmd = $update_genome_script." \"$name\" ";
    print "$species_specific_cmd\n";
    system("$species_specific_cmd") if $write;
    $all_species{$name} = 1;
    $species_count++ if !exists $all_species{$name};;
#    $species_count_hash{'EnsEMBL'}{$name} = 1;
}
}
########################################################
# Load non-Ensembl species into genome_db
########################################################
my $insert_command = "mysql -h $db_host -P $db_port -u $db_user -p$db_pass -D $master_db_name -e 'INSERT INTO genome_db (taxon_id, name, assembly) VALUES ";
# read species file
open my $species_file, "<$species_file" or die "\tcould not read species file $species_file\n";
while(<$species_file>){
	chomp;
	next if /^$/;
    my ($aa_file,$cds_file,$ncbi_tax_id,$name, $assembly) = split(/\s+/, $_);
	my $specific_insert_command = $insert_command." ($ncbi_tax_id, \"$name\", \"$assembly\") '";
	print "$specific_insert_command\n";
    system("$specific_insert_command") if $write;
    $all_species{$name} = 1;
    $species_count++ if !exists $all_species{$name};;
}

print "Added ".keys(%all_species)." species so far\n";

foreach my $species(keys(%all_species)){
    print "SPECIES: $species\n";
}
# For species from ensembl metazoa
# perl backup/update_genome.pl --reg_conf $ENSEMBL_CVS_ROOT_DIR/treefam/scripts/pipeline/production_treefam_reg_conf.pl --compara compara_master --species "Amphimedon queenslandica" --offset 100000000  --taxon_id 400682perl backup/update_genome.pl --reg_conf $ENSEMBL_CVS_ROOT_DIR/treefam/scripts/pipeline/production_treefam_reg_conf.pl --compara compara_master --species "Trichoplax adhaerens" --offset 10000000  --taxon_id 10228# Add Genome_db entries for fasta-species
# INSERT INTO genome_db (taxon_id, name) VALUES (81824, "monosiga_brevicollis");
# INSERT INTO genome_db (taxon_id, name) VALUES (6412, "helobdella_robusta");
# INSERT INTO genome_db (taxon_id, name) VALUES (73382, "capitella_teleta");
# INSERT INTO genome_db (taxon_id, name) VALUES (218847, "proterospongia");
# # Now correctly populate mlss tables
# # Need to store all genome_db_ids in a variable
# 
 my $export_all_genomes_cmd = "export ALL_GENOMEDB_IDS=`mysql -h $db_host -P $db_port -u $db_user -p$db_pass $master_db_name -N -e \"select group_concat(genome_db_id order by genome_db_id) from genome_db where (assembly_default=1 and taxon_id is not NULL)\" | cat`";
 print $export_all_genomes_cmd."\n";
system($export_all_genomes_cmd) if $write_orthologs;
# ## choose a temp. directory where the output will be generated:
# 
my $export_dir = 'export MLSS_DIR="/tmp/mlss_creation"';
print $export_dir."\n";
system($export_dir) if $write_orthologs;
print "mkdir \$MLSS_DIR\n";
mkdir("/tmp/mlss_creation") if $write_orthologs;
#     export MLSS_DIR="/tmp/mlss_creation"
#     mkdir $MLSS_DIR
# 
# # Make sure that .profile is read
# ## run the loading script several times:
# 
# # orthologues
 my $orthologs_cmd = " echo -e \"201\\n\" | perl \$ENSEMBL_CVS_ROOT_DIR/ensembl-compara/scripts/pipeline/create_mlss.pl --f --reg_conf $registry_file --pw --genome_db_id \"\$ALL_GENOMEDB_IDS\" 1>\$MLSS_DIR/create_mlss.ENSEMBL_ORTHOLOGUES.201.out 2>\$MLSS_DIR/create_mlss.ENSEMBL_ORTHOLOGUES.201.err";
 print $orthologs_cmd."\n";
system($orthologs_cmd) if $write_orthologs;
# 
# # paralogues btw
my $paralogs_btw_cmd = " echo -e \"202\\n\" | perl \$ENSEMBL_CVS_ROOT_DIR/ensembl-compara/scripts/pipeline/create_mlss.pl --f \
 --reg_conf $registry_file \
 --pw --genome_db_id \"\$ALL_GENOMEDB_IDS\" 1>\$MLSS_DIR/create_mlss.ENSEMBL_PARALOGUES.btw.202.out 2>\$MLSS_DIR/create_mlss.ENSEMBL_PARALOGUES.btw.202.err";
print $paralogs_btw_cmd."\n";
system($paralogs_btw_cmd) if $write_orthologs;
 
# # paralogues wth
my $paralogs_wth_cmd =  "echo -e \"202\\n\" | perl \$ENSEMBL_CVS_ROOT_DIR/ensembl-compara/scripts/pipeline/create_mlss.pl --f \
 --reg_conf $registry_file \
 --sg --genome_db_id \"\$ALL_GENOMEDB_IDS\" 1>\$MLSS_DIR/create_mlss.ENSEMBL_PARALOGUES.wth.202.out 2>\$MLSS_DIR/create_mlss.ENSEMBL_PARALOGUES.wth.202.err";
print "$paralogs_wth_cmd\n";
system($paralogs_wth_cmd) if $write_orthologs;
 # 
# # proteintrees
my $protein_trees_cmd = "echo -e \"401\\n\" | perl \$ENSEMBL_CVS_ROOT_DIR/ensembl-compara/scripts/pipeline/create_mlss.pl --f \
 --reg_conf $registry_file \
 --name \"protein trees\" --genome_db_id \"\$ALL_GENOMEDB_IDS\" 1>\$MLSS_DIR/create_mlss.PROTEIN_TREES.401.out 2>\$MLSS_DIR/create_mlss.PROTEIN_TREES.401.err";

 print "$protein_trees_cmd\n";
system($protein_trees_cmd) if $write_orthologs;
 #
### Now, make sure the pipeline can find the non-ensembl species
# Check config-file for location of json-file "/lustre/scratch109/sanger/fs9/treefam/species_list.json"
## make appropriate changes to production script
#
#

# make sure the HMMS are in the right format (hmmer 2) and that there is a consensus file
# script in "/nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/treefam_tools/make_hmms_ready/run_all.sh"


#
