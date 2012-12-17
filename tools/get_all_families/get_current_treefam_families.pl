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


my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(-user   => 'treefam_user',
                                 #-dbname => 'treefam_homology_67hmm',
                                 -dbname => 'treefam_production_9_68',
                                                                -host   => 'web-mei-treefam',
                                                                -pass => 'readonly',
                                                                -port => '3365');
#my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
#    -user   => 'anonymous',
#    -dbname => 'ensembl_compara_67',
#    -host   => 'ensembldb.ensembl.org',
#    -port   => '5306'
#);

print "read all treefam families\n";
my $families_file = "familyA.txt.table";
my @all_families = `cat $families_file`;
my %treefam8_hash;
foreach(@all_families){
    chomp;
    my ($fam,$symbol,$desc,$prevS,$prevF) = split(/\t/);
    $treefam8_hash{$fam}{'symbol'} = $symbol;
    $desc =~ s/'|"//g;
    $treefam8_hash{$fam}{'desc'} = $desc;
    $treefam8_hash{$fam}{'prevS'} = $prevS;
    $treefam8_hash{$fam}{'prevF'} = $prevF;
}

my $genetree_adaptor = $db->get_GeneTreeAdaptor;
my $tree = $genetree_adaptor->fetch_by_stable_id("TF101003");
my @all_trees = @{$genetree_adaptor->fetch_all( -CLUSTERSET_ID => "default")};
#my @all_trees; 
push(@all_trees,$tree);
print "Found ".scalar(@all_trees)." trees\n";
my $familyCount = 0;
my $family_file = "treefam9_families.json";
open my $file_out,">", $family_file or die "Could not open file\n";
print {$file_out} "[\n";
foreach my $gt(@all_trees){
    my %gt_names;
    my $tagvalue_hashref = $gt->get_tagvalue_hash();
    if(!keys(%$tagvalue_hashref)){
        die "Could not get tagvalue_hashref for tree\n";
    }
	my ($modelName,$alnPercentIdentity,$alnLength,$geneCount, $pfam) = (
		$tagvalue_hashref->{model_name},
		$tagvalue_hashref->{aln_percent_identity},
		$tagvalue_hashref->{aln_length},
		$tagvalue_hashref->{gene_count},
		$tagvalue_hashref->{pfam},

	);
    ## get taxonomic distribution
    ## pfam domains?
    next if $modelName !~ /TF1/;
    #next if $modelName !~ /TF101003/;
	$alnPercentIdentity = ($alnPercentIdentity eq "")? "NaN": int($alnPercentIdentity);
	$alnLength = ($alnLength eq "")? "NaN": $alnLength;
	$geneCount = ($geneCount eq "")? "NaN": $geneCount;
    my @pfams;
    if(defined($pfam) && $pfam ne ''){
         my @array = split(" ",$pfam);         
         foreach(@array){
            my @arr2 = split("->",$_);
            push(@pfams, $arr2[0]);
        }
    }
    my $pfams = (scalar(@pfams))? join(",",@pfams): "No hits";
    $gt_names{$modelName}{geneCount} = ($geneCount eq "")? "NaN": $geneCount;
	
    #print "Checking $modelName\n";
    #my $taxa_counts = &get_taxa_counts(($gt->root)->nhx_format;);
    my $gt_root = $gt->root;
    
    my $root_taxon = "NaN";
    if($gt_root->has_tag('taxon_name')){
        $root_taxon = $gt_root->get_value_for_tag('taxon_name');
    }

    print {$file_out}  "{\"modelName\":\"$modelName\", \"hgncSymbol\":\"".$treefam8_hash{$modelName}{symbol}."\", \"geneCount\":\"$geneCount\", \"percentIdentity\":\"$alnPercentIdentity\",\"alnLength\":\"$alnLength\",\"rootTaxon\":\"$root_taxon\",\"description\":\"".$treefam8_hash{$modelName}{desc}."\"},\n";
	#last if $familyCount++ > 300;
}

print {$file_out} "]";
close $file_out;

exit;

