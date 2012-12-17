#!/usr/bin/env perl -w
use LWP::Simple;
use File::Basename;
use Data::Dumper;
my $usage = "Usage: $0 source\n";

my $source = "ensembl";
my $source_transcriptDir = '/lustre/scratch109/sanger/fs9/treefam/all_ensembl_sequences/transcripts';
my $source_proteinDir = '/lustre/scratch109/sanger/fs9/treefam/all_ensembl_sequences/proteins';
my $source_gtfDir = '/lustre/scratch109/sanger/fs9/treefam/all_ensembl_sequences/gtf';
my $release = '69';


mkdir($source_transcriptDir ) if !-e $source_transcriptDir;
mkdir($source_proteinDir) if !-e $source_proteinDir;
mkdir($source_gtfDir) if !-e $source_gtfDir;
my $warehouse_source_dir = "/warehouse/pfam01/treefam/working_dir/sequences/";

my $current_ensemblDir = "ftp://ftp.ensembl.org/pub/release-$release/";
my $species_list = "species2.list";
my $cnda_file_ending = "cdna.all.fa.gz";
my $protein_file_ending = "pep.all.fa.gz";
my $gtf_file_ending = "gtf.gz";

my $local_file_ending = "fa.gz";


#### download all




if($source eq 'ensembl'){
my %species2shortcut = (
"ailuropoda_melanoleuca" => {'name' => 'Ailuropoda_melanoleuca','shortcut' => 'ailMel1', 'swcode' => 'AILME', 'taxid' => '9646'},
"anolis_carolinensis" => {'name' => 'Anolis_carolinensis','shortcut' => 'AnoCar2.0', 'swcode' => 'ANOCA', 'taxid' => '28377'},
"bos_taurus" => {'name' => 'Bos_taurus','shortcut' => 'UMD3.1', 'swcode' => 'BOVIN', 'taxid' => '9913'},
#"caenorhabditis_elegans" => {'name' => 'Caenorhabditis_elegans','shortcut' => 'WS220', 'swcode' => 'CAEEL', 'taxid' => '6239'},
"callithrix_jacchus" => {'name' => 'Callithrix_jacchus','shortcut' => 'C_jacchus3.2.1', 'swcode' => 'CALJA', 'taxid' => '9483'},
"canis_familiaris" => {'name' => 'Canis_familiaris', 'shortcut' => 'CanFam3.1', 'swcode' => 'CANFA', 'taxid' => '9615'},
"cavia_porcellus" => {'name' => 'Cavia_porcellus', 'shortcut' => 'cavPor3', 'swcode' => 'CAVPO', 'taxid' => '10141'},
"choloepus_hoffmanni" => {'name' => 'Choloepus_hoffmanni', 'shortcut' => 'choHof1', 'swcode' => 'CHOHO', 'taxid' => '9358'},
"ciona_intestinalis" => {'name' => 'Ciona_intestinalis', 'shortcut' => 'KH', 'swcode' => 'CIOIN', 'taxid' => '7719'},
"ciona_savignyi" => {'name' => 'Ciona_savignyi', 'shortcut' => 'CSAV2.0', 'swcode' => 'CIOSA', 'taxid' => '51511'},
"danio_rerio" => {'name' => 'Danio_rerio', 'shortcut' => 'Zv9', 'swcode' => 'BRARE', 'taxid' => '7955'},
"dasypus_novemcinctus" => {'name' => 'Dasypus_novemcinctus', 'shortcut' => 'dasNov2', 'swcode' => 'DASNO', 'taxid' => '9361'},
"dipodomys_ordii" => {'name' => 'Dipodomys_ordii', 'shortcut' => 'dipOrd1', 'swcode' => 'DIPOR', 'taxid' => '10020'},
#"drosophila_melanogaster" => {'name' => 'Drosophila_melanogaster', 'shortcut' => 'BDGP5', 'swcode' => 'DROME', 'taxid' => '7227'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP5.66.cdna.all.fa.gz
"echinops_telfairi" => {'name' => 'Echinops_telfairi', 'shortcut' => 'TENREC', 'swcode' => 'ECHTE', 'taxid' => '9371'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/echinops_telfairi/cdna/Echinops_telfairi.TENREC.66.cdna.all.fa.gz
"equus_caballus" => {'name' => 'Equus_caballus', 'shortcut' => 'EquCab2', 'swcode' => 'HORSE', 'taxid' => '9796'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/equus_caballus/cdna/Equus_caballus.EquCab2.66.cdna.all.fa.gz
"erinaceus_europaeus" => {'name' => 'Erinaceus_europaeus', 'shortcut' => 'HEDGEHOG', 'swcode' => 'ERIEU', 'taxid' => '9365'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/erinaceus_europaeus/cdna/Erinaceus_europaeus.HEDGEHOG.66.cdna.all.fa.gz
"felis_catus" => {'name' => 'Felis_catus', 'shortcut' => 'CAT', 'swcode' => 'FELCA', 'taxid' => '9685'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/felis_catus/cdna/Felis_catus.CAT.66.cdna.all.fa.gz
"gadus_morhua" => {'name' => 'Gadus_morhua', 'shortcut' => 'gadMor1', 'swcode' => 'GADMO', 'taxid' => '8049'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/gadus_morhua/cdna/Gadus_morhua.gadMor1.66.cdna.all.fa.gz
"gallus_gallus" => {'name' => 'Gallus_gallus', 'shortcut' => 'WASHUC2', 'swcode' => 'CHICK', 'taxid' => '9031'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/gallus_gallus/cdna/Gallus_gallus.WASHUC2.66.cdna.all.fa.gz
"gasterosteus_aculeatus" => {'name' => 'Gasterosteus_aculeatus', 'shortcut' => 'BROADS1', 'swcode' => 'GASAC', 'taxid' => '69293'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/gasterosteus_aculeatus/cdna/Gasterosteus_aculeatus.BROADS1.66.cdna.all.fa.gz
"gorilla_gorilla" => {'name' => 'Gorilla_gorilla', 'shortcut' => 'gorGor3.1', 'swcode' => 'GORGO', 'taxid' => '9595'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/gorilla_gorilla/cdna/Gorilla_gorilla.gorGor3.1.66.cdna.all.fa.gz
"homo_sapiens" => {'name' => 'Homo_sapiens', 'shortcut' => 'GRCh37', 'swcode' => 'HUMAN', 'taxid' => '9606'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.66.cdna.all.fa.gz
"ictidomys_tridecemlineatus" => {'name' => 'Ictidomys_tridecemlineatus', 'shortcut' => 'spetri2', 'swcode' => 'LATCH', 'taxid' => '7897'},
"latimeria_chalumnae" => {'name' => 'Latimeria_chalumnae', 'shortcut' => 'LatCha1', 'swcode' => 'LATCH', 'taxid' => '7897'},
"loxodonta_africana" => {'name' => 'Loxodonta_africana', 'shortcut' => 'loxAfr3', 'swcode' => 'LOXAF', 'taxid' => '9785'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/loxodonta_africana/cdna/Loxodonta_africana.loxAfr3.66.cdna.all.fa.gz
"macaca_mulatta" => {'name' => 'Macaca_mulatta', 'shortcut' => 'MMUL_1', 'swcode' => 'MACMU', 'taxid' => '9544'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/macaca_mulatta/cdna/Macaca_mulatta.MMUL_1.66.cdna.all.fa.gz
"macropus_eugenii" => {'name' => 'Macropus_eugenii', 'shortcut' => 'Meug_1.0', 'swcode' => 'MACEU', 'taxid' => '9315'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/macropus_eugenii/cdna/Macropus_eugenii.Meug_1.0.66.cdna.all.fa.gz
"meleagris_gallopavo" => {'name' => 'Meleagris_gallopavo', 'shortcut' => 'UMD2', 'swcode' => 'MELGA', 'taxid' => '9103'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/meleagris_gallopavo/cdna/Meleagris_gallopavo.UMD2.66.cdna.all.fa.gz
"microcebus_murinus" => {'name' => 'Microcebus_murinus', 'shortcut' => 'micMur1', 'swcode' => 'MICMU', 'taxid' => '30608'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/microcebus_murinus/cdna/Microcebus_murinus.micMur1.66.cdna.all.fa.gz
"monodelphis_domestica" => {'name' => 'Monodelphis_domestica', 'shortcut' => 'BROADO5', 'swcode' => 'MONDO', 'taxid' => '13616'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/monodelphis_domestica/cdna/Monodelphis_domestica.BROADO5.66.cdna.all.fa.gz
"mus_musculus" => {'name' => 'Mus_musculus', 'shortcut' => 'GRCm38', 'swcode' => 'MOUSE', 'taxid' => '10090'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/mus_musculus/cdna/Mus_musculus.NCBIM37.66.cdna.all.fa.gz
"myotis_lucifugus" => {'name' => 'Myotis_lucifugus', 'shortcut' => 'Myoluc2.0', 'swcode' => 'MYOLU', 'taxid' => '59463'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/myotis_lucifugus/cdna/Myotis_lucifugus.Myoluc2.0.66.cdna.all.fa.gz
"nomascus_leucogenys" => {'name' => 'Nomascus_leucogenys', 'shortcut' => 'Nleu1.0', 'swcode' => 'NOLEU', 'taxid' => '61853'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/nomascus_leucogenys/cdna/Nomascus_leucogenys.Nleu1.0.66.cdna.all.fa.gz
"ochotona_princeps" => {'name' => 'Ochotona_princeps', 'shortcut' => 'pika', 'swcode' => 'OCHPR', 'taxid' => 'OCHPR'},
"oreochromis_niloticus" => {'name' => 'Oreochromis_niloticus', 'shortcut' => 'Orenil1.0', 'swcode' => 'OCHPR', 'taxid' => '8128'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/ochotona_princeps/cdna/Ochotona_princeps.pika.66.cdna.all.fa.gz
"ornithorhynchus_anatinus" => {'name' => 'Ornithorhynchus_anatinus', 'shortcut' => 'OANA5', 'swcode' => 'ORNAN', 'taxid' => '9258'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/ornithorhynchus_anatinus/cdna/Ornithorhynchus_anatinus.OANA5.66.cdna.all.fa.gz
"oryctolagus_cuniculus" => {'name' => 'Oryctolagus_cuniculus', 'shortcut' => 'oryCun2', 'swcode' => 'RABIT', 'taxid' => '9986'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/oryctolagus_cuniculus/cdna/Oryctolagus_cuniculus.oryCun2.66.cdna.all.fa.gz
"oryzias_latipes" => {'name' => 'Oryzias_latipes', 'shortcut' => 'MEDAKA1', 'swcode' => 'ORYLA', 'taxid' => '8090'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/oryzias_latipes/cdna/Oryzias_latipes.MEDAKA1.66.cdna.all.fa.gz
"otolemur_garnettii" => {'name' => 'Otolemur_garnettii', 'shortcut' => 'OtoGar3', 'swcode' => 'OTOGA', 'taxid' => '30611'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/otolemur_garnettii/cdna/Otolemur_garnettii.OtoGar3.66.cdna.all.fa.gz
"pan_troglodytes" => {'name' => 'Pan_troglodytes', 'shortcut' => 'CHIMP2.1.4', 'swcode' => 'PANTR', 'taxid' => '51511'},

"pelodiscus_sinensis" => {'name' => 'Pelodiscus_sinensis', 'shortcut' => 'PelSin_1.0', 'swcode' => 'PANTR', 'taxid' => '13735'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/pan_troglodytes/cdna/Pan_troglodytes.CHIMP2.1.4.66.cdna.all.fa.gz
"petromyzon_marinus" => {'name' => 'Petromyzon_marinus', 'shortcut' => 'Pmarinus_7.0', 'swcode' => 'PETMA', 'taxid' => '7757'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/petromyzon_marinus/cdna/Petromyzon_marinus.Pmarinus_7.0.66.cdna.all.fa.gz
"pongo_abelii" => {'name' => 'Pongo_abelii', 'shortcut' => 'PPYG2', 'swcode' => 'PONAM', 'taxid' => '62124'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/pongo_abelii/cdna/Pongo_abelii.PPYG2.66.cdna.all.fa.gz
"procavia_capensis" => {'name' => 'Procavia_capensis', 'shortcut' => 'proCap1', 'swcode' => 'PROCA', 'taxid' => '9813'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/procavia_capensis/cdna/Procavia_capensis.proCap1.66.cdna.all.fa.gz
"pteropus_vampyrus" => {'name' => 'Pteropus_vampyrus', 'shortcut' => 'pteVam1', 'swcode' => 'PTEVA', 'taxid' => '132908'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/pteropus_vampyrus/cdna/Pteropus_vampyrus.pteVam1.66.cdna.all.fa.gz
"rattus_norvegicus" => {'name' => 'Rattus_norvegicus', 'shortcut' => 'RGSC3.4', 'swcode' => 'RAT', 'taxid' => '10116'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.RGSC3.4.66.cdna.all.fa.gz
"saccharomyces_cerevisiae" => {'name' => 'Saccharomyces_cerevisiae', 'shortcut' => 'EF4', 'swcode' => 'YEAST', 'taxid' => '4932'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.EF4.66.cdna.all.fa.gz
"sarcophilus_harrisii" => {'name' => 'Sarcophilus_harrisii', 'shortcut' => 'DEVIL7.0', 'swcode' => 'SARHA', 'taxid' => '9305'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/sarcophilus_harrisii/cdna/Sarcophilus_harrisii.DEVIL7.0.66.cdna.all.fa.gz
"sorex_araneus" => {'name' => 'Sorex_araneus', 'shortcut' => 'COMMON_SHREW1', 'swcode' => 'SORAR', 'taxid' => '42254'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/sorex_araneus/cdna/Sorex_araneus.COMMON_SHREW1.66.cdna.all.fa.gz
#"spermophilus_tridecemlineatus" => {'name' => 'Spermophilus_tridecemlineatus', 'shortcut' => 'SQUIRREL', 'swcode' => 'SPETR', 'taxid' => '43179'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/spermophilus_tridecemlineatus/cdna/Spermophilus_tridecemlineatus.SQUIRREL.66.cdna.all.fa.gz
"sus_scrofa" => {'name' => 'Sus_scrofa', 'shortcut' => 'Sscrofa10.2', 'swcode' => 'PIG', 'taxid' => '9823'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/sus_scrofa/cdna/Sus_scrofa.Sscrofa9.66.cdna.all.fa.gz
"taeniopygia_guttata" => {'name' => 'Taeniopygia_guttata', 'shortcut' => 'taeGut3.2.4', 'swcode' => 'TAGUT', 'taxid' => '59729'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/taeniopygia_guttata/cdna/Taeniopygia_guttata.taeGut3.2.4.66.cdna.all.fa.gz
"takifugu_rubripes" => {'name' => 'Takifugu_rubripes', 'shortcut' => 'FUGU4', 'swcode' => 'FUGU', 'taxid' => '31033'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/takifugu_rubripes/cdna/Takifugu_rubripes.FUGU4.66.cdna.all.fa.gz
"tarsius_syrichta" => {'name' => 'Tarsius_syrichta', 'shortcut' => 'tarSyr1', 'swcode' => 'TARSY', 'taxid' => '9478'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/tarsius_syrichta/cdna/Tarsius_syrichta.tarSyr1.66.cdna.all.fa.gz
"tetraodon_nigroviridis" => {'name' => 'Tetraodon_nigroviridis', 'shortcut' => 'TETRAODON8', 'swcode' => 'TETNG', 'taxid' => '99883'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/tetraodon_nigroviridis/cdna/Tetraodon_nigroviridis.TETRAODON8.66.cdna.all.fa.gz
"tupaia_belangeri" => {'name' => 'Tupaia_belangeri', 'shortcut' => 'TREESHREW', 'swcode' => 'TUPGB', 'taxid' => '37347'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/tupaia_belangeri/cdna/Tupaia_belangeri.TREESHREW.66.cdna.all.fa.gz
"tursiops_truncatus" => {'name' => 'Tursiops_truncatus', 'shortcut' => 'turTru1', 'swcode' => 'TURTR', 'taxid' => '9739'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/tursiops_truncatus/cdna/Tursiops_truncatus.turTru1.66.cdna.all.fa.gz
"vicugna_pacos" => {'name' => 'Vicugna_pacos', 'shortcut' => 'vicPac1', 'swcode' => 'VIPAC', 'taxid' => '30538'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/vicugna_pacos/cdna/Vicugna_pacos.vicPac1.66.cdna.all.fa.gz
"xenopus_tropicalis" => {'name' => 'Xenopus_tropicalis', 'shortcut' => 'JGI_4.2', 'swcode' => 'XENTR', 'taxid' => '8364'},
#ftp://ftp.ensembl.org/pub/release-66/fasta/xenopus_tropicalis/cdna/Xenopus_tropicalis.JGI_4.2.66.cdna.all.fa.gz
);

print "\thave ".keys(%species2shortcut)." species\n";
	$current_ensemblDir = "ftp://ftp.ensembl.org/pub/release-$release/";
	my $speciesListString;
	foreach my $species_name(keys(%species2shortcut))
	{
		#'cdna' => 'ftp://ftp.ensembl.org/pub/release-66/fasta/ailuropoda_melanoleuca/cdna/Ailuropoda_melanoleuca.ailMel1.66.cdna.all.fa.gz',
		#'protein' => 'ftp://ftp.ensembl.org/pub/release-66/fasta/ailuropoda_melanoleuca/pep/Ailuropoda_melanoleuca.ailMel1.66.pep.all.fa.gz',
		#'gtf' => 'ftp://ftp.ensembl.org/pub/release-66/gtf/ailuropoda_melanoleuca/Ailuropoda_melanoleuca.ailMel1.66.gtf.gz',
		print "\tDownloading $species_name\n";
		#my $cndaSource = $species{$species_name}{'cdna'};
		my $speciesName2 = $species2shortcut{$species_name}{'name'};
		my $shortcut = $species2shortcut{$species_name}{'shortcut'};
		my $swcode = $species2shortcut{$species_name}{'swcode'};
		my $taxid = $species2shortcut{$species_name}{'taxid'};
		
		my $fileName = $speciesName2.".".$shortcut.".".$release.".".$cnda_file_ending;
		my $cndaSource = $current_ensemblDir."/fasta/".$species_name."/cdna/".$fileName; # remote
		my $cndaFile = $source_transcriptDir."/".$fileName; # local 

		unless(-e $cndaFile)
		{
			print "\tDownloading $cndaSource ($cndaFile)\n";
			getstore($cndaSource, $cndaFile) or die "CHECK:$cndaSource\n";
		}	
		#my $proteinSource = $species{$species_name}{'protein'};
		$fileName = $speciesName2.".".$shortcut.".".$release.".".$protein_file_ending;
		my $proteinSource = $current_ensemblDir."/fasta/".$species_name."/pep/".$fileName; #remote
		#my $proteinFile = $source_proteinDir."/".$fileName; # local 
		my $proteinFile = $source_proteinDir."/".$speciesName2.".".$local_file_ending  ; # local 
		unless(-e $proteinFile)
		{
			print "\tDownloading $proteinSource\n";
			getstore($proteinSource, $proteinFile) or die "CHECK:$proteinSource\n";
		}
		$fileName =~ s/\.gz//;
		$speciesName2toprint = $speciesName2;
		$speciesName2toprint =~ s/_/ /;
		$speciesListString .= "$speciesName2toprint\t$fileName\t$swcode\t$taxid\tENS\n";
		
		$fileName = $speciesName2.".".$shortcut.".".$release.".".$gtf_file_ending;
		#my $gtfFile = $source_gtfDir."/".basename($species{$species_name}{'gtf'});
		my $gtfFile = $source_gtfDir."/".$fileName;
		my $gtfSource = $current_ensemblDir."/gtf/".$species_name."/".$fileName;
		
		unless(-e $gtfFile){
			print "\tDownloading $gtfSource\n";
			getstore($gtfSource, $gtfFile) or die "CHECK:$gtfSource\n";
		}	
	}
	
	## now write species list
	# Ailuropoda melanoleuca	Ailuropoda_melanoleuca.ailMel1.66.pep.all.fa	AILME	9646	ENS	
	# Acropora digitifera	adi_v1.0.1.prot.fa	ACRDI	70779	OTH
	
	write_to_file({file_name => $species_list, text => $speciesListString});
}
if($source eq 'ensembl-metazoa'){
	## metazoa uses different releases
	$release = '15';
	
	my %species2shortcut = (
	"acyrthosiphon_pisum" => {'name' => 'Acyrthosiphon_pisum','shortcut' => 'Acyr_2.0', 'swcode' => 'ACYPI', 'taxid' => '7029'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/acyrthosiphon_pisum/cdna/Acyrthosiphon_pisum.Acyr2.12.cdna.all.fa.gz
"aedes_aegypti" => {'name' => 'Aedes_aegypti','shortcut' => 'AaegL1', 'swcode' => 'AEDAE', 'taxid' => '7159'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/aedes_aegypti/cdna/Aedes_aegypti.AaegL1.12.cdna.all.fa.gz
"amphimedon_queenslandica" => {'name' => 'Amphimedon_queenslandica','shortcut' => 'Aqu1', 'swcode' => 'AMQUE', 'taxid' => '400682'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/amphimedon_queenslandica/cdna/Amphimedon_queenslandica.Aqu1.12.cdna.all.fa.gz
"anopheles_gambiae" => {'name' => 'Anopheles_gambiae','shortcut' => 'AgamP3', 'swcode' => 'ANOGA', 'taxid' => '7165'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/anopheles_gambiae/cdna/Anopheles_gambiae.AgamP3.12.cdna.all.fa.gz
"apis_mellifera" => {'name' => 'Apis_mellifera','shortcut' => 'Amel_2.0', 'swcode' => 'APIME', 'taxid' => '7460'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/apis_mellifera/cdna/Apis_mellifera.Amel_2.0.12.cdna.all.fa.gz
"atta_cephalotes" => {'name' => 'Atta_cephalotes','shortcut' => 'Attacep1.0', 'swcode' => 'ATCEP', 'taxid' => '12957'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/atta_cephalotes/cdna/Atta_cephalotes.Acep1.12.cdna.all.fa.gz
"bombyx_mori" => {'name' => 'Bombyx_mori','shortcut' => 'Bmor1', 'swcode' => 'BOMMO', 'taxid' => '7091'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/bombyx_mori/cdna/Bombyx_mori.Bmor1.12.cdna.all.fa.gz
"caenorhabditis_brenneri" => {'name' => 'Caenorhabditis_brenneri','shortcut' => 'CB601', 'swcode' => 'CABRE', 'taxid' => '135651'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/caenorhabditis_brenneri/cdna/Caenorhabditis_brenneri.CB601.12.cdna.all.fa.gz
"caenorhabditis_briggsae" => {'name' => 'Caenorhabditis_briggsae','shortcut' => 'CB4', 'swcode' => 'CAEBR', 'taxid' => '6238'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/caenorhabditis_briggsae/cdna/Caenorhabditis_briggsae.CB3.12.cdna.all.fa.gz
"caenorhabditis_elegans" => {'name' => 'Caenorhabditis_elegans','shortcut' => 'WBcel215', 'swcode' => 'CAEEL', 'taxid' => '6239'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WS220.12.cdna.all.fa.gz
"caenorhabditis_japonica" => {'name' => 'Caenorhabditis_japonica','shortcut' => 'CJ302', 'swcode' => 'CAEJA', 'taxid' => '281687'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/caenorhabditis_japonica/cdna/Caenorhabditis_japonica.CJ302.12.cdna.all.fa.gz
"caenorhabditis_remanei" => {'name' => 'Caenorhabditis_remanei','shortcut' => 'CR2', 'swcode' => 'CAERE', 'taxid' => '31234'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/caenorhabditis_remanei/cdna/Caenorhabditis_remanei.CR2.12.cdna.all.fa.gz
"culex_quinquefasciatus" => {'name' => 'Culex_quinquefasciatus','shortcut' => 'CpipJ1', 'swcode' => 'CULQU', 'taxid' => '7176'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/culex_quinquefasciatus/cdna/Culex_quinquefasciatus.CpipJ1.12.cdna.all.fa.gz
"danaus_plexippus" => {'name' => 'Danaus_plexippus','shortcut' => 'DanPle_1.0', 'swcode' => 'DAPPU', 'taxid' => '13037'},
"daphnia_pulex" => {'name' => 'Daphnia_pulex','shortcut' => 'Dappu1', 'swcode' => 'DAPPU', 'taxid' => '6669'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/daphnia_pulex/cdna/Daphnia_pulex.Dappu1.12.cdna.all.fa.gz
"drosophila_ananassae" => {'name' => 'Drosophila_ananassae','shortcut' => 'dana_caf1', 'swcode' => 'DROAN', 'taxid' => '7217'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/drosophila_ananassae/cdna/Drosophila_ananassae.dana_r1.3_FB2008_07.12.cdna.all.fa.gz
"drosophila_erecta" => {'name' => 'Drosophila_erecta','shortcut' => 'dere_caf1', 'swcode' => 'DROER', 'taxid' => '7220'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/drosophila_erecta/cdna/Drosophila_erecta.dere_r1.3_FB2008_07.12.cdna.all.fa.gz
"drosophila_grimshawi" => {'name' => 'Drosophila_grimshawi','shortcut' => 'dgri_caf1', 'swcode' => 'DROGR', 'taxid' => '7222'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/drosophila_grimshawi/cdna/Drosophila_grimshawi.dgri_r1.3_FB2008_07.12.cdna.all.fa.gz
"drosophila_melanogaster" => {'name' => 'Drosophila_melanogaster','shortcut' => 'BDGP5', 'swcode' => 'DROME', 'taxid' => '7227'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP5.12.cdna.all.fa.gz
"drosophila_mojavensis" => {'name' => 'Drosophila_mojavensis','shortcut' => 'dmoj_caf1', 'swcode' => 'DROMO', 'taxid' => '7230'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/drosophila_mojavensis/cdna/Drosophila_mojavensis.dmoj_r1.3_FB2008_07.12.cdna.all.fa.gz
"drosophila_persimilis" => {'name' => 'Drosophila_persimilis','shortcut' => 'dper_caf1', 'swcode' => 'DROPE', 'taxid' => '7234'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/drosophila_persimilis/cdna/Drosophila_persimilis.dper_r1.3_FB2008_07.12.cdna.all.fa.gz
"drosophila_pseudoobscura" => {'name' => 'Drosophila_pseudoobscura','shortcut' => 'HGSC2', 'swcode' => 'DROPS', 'taxid' => '7237'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/drosophila_pseudoobscura/cdna/Drosophila_pseudoobscura.HGSC2.12.cdna.all.fa.gz
"drosophila_sechellia" => {'name' => 'Drosophila_sechellia','shortcut' => 'dsec_caf1', 'swcode' => 'DROSE', 'taxid' => '7238'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/drosophila_sechellia/cdna/Drosophila_sechellia.dsec_r1.3_FB2008_07.12.cdna.all.fa.gz
"drosophila_simulans" => {'name' => 'Drosophila_simulans','shortcut' => 'dsim_r1.3_FB2008_07', 'swcode' => 'DROSI', 'taxid' => '7240'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/drosophila_simulans/cdna/Drosophila_simulans.dsim_r1.3_FB2008_07.12.cdna.all.fa.gz
"drosophila_virilis" => {'name' => 'Drosophila_virilis','shortcut' => 'dvir_caf1', 'swcode' => 'DROVI', 'taxid' => '7244'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/drosophila_virilis/cdna/Drosophila_virilis.dvir_r1.2_FB2008_07.12.cdna.all.fa.gz
"drosophila_willistoni" => {'name' => 'Drosophila_willistoni','shortcut' => 'dwil_caf1', 'swcode' => 'DROWI', 'taxid' => '7260'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/drosophila_willistoni/cdna/Drosophila_willistoni.dwil_r1.3_FB2008_07.12.cdna.all.fa.gz
"drosophila_yakuba" => {'name' => 'Drosophila_yakuba','shortcut' => 'dyak_r1.3_FB2008_07', 'swcode' => 'DROYA', 'taxid' => '7245'},
"heliconius_melpomene" => {'name' => 'Heliconius_melpomene','shortcut' => 'Hmel1', 'swcode' => 'IXOSC', 'taxid' => '34740'},
#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-15/fasta/heliconius_melpomene/pep/Heliconius_melpomene.Hmel1.15.pep.all.fa.gz
"ixodes_scapularis" => {'name' => 'Ixodes_scapularis','shortcut' => 'IscaW1', 'swcode' => 'IXOSC', 'taxid' => '6945'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/ixodes_scapularis/cdna/Ixodes_scapularis.IscaW1.12.cdna.all.fa.gz
"nematostella_vectensis" => {'name' => 'Nematostella_vectensis','shortcut' => 'ASM20922v1', 'swcode' => 'NEMVE', 'taxid' => '45351'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/nematostella_vectensis/cdna/Nematostella_vectensis.Nemve1.12.cdna.all.fa.gz
"pediculus_humanus" => {'name' => 'Pediculus_humanus','shortcut' => 'PhumU1', 'swcode' => 'PEDHU', 'taxid' => '121225'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/pediculus_humanus/cdna/Pediculus_humanus.PhumU1.12.cdna.all.fa.gz
"pristionchus_pacificus" => {'name' => 'Pristionchus_pacificus','shortcut' => 'pp1', 'swcode' => 'PRIPA', 'taxid' => '54126'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/pristionchus_pacificus/cdna/Pristionchus_pacificus.pp1.12.cdna.all.fa.gz
"schistosoma_mansoni" => {'name' => 'Schistosoma_mansoni','shortcut' => 'sma_v3.1', 'swcode' => 'SCHMA', 'taxid' => '6183'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/schistosoma_mansoni/cdna/Schistosoma_mansoni.sma_v3.1.12.cdna.all.fa.gz
"strongylocentrotus_purpuratus" => {'name' => 'Strongylocentrotus_purpuratus','shortcut' => 'Spur2.5', 'swcode' => 'STRPU', 'taxid' => '7668'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/strongylocentrotus_purpuratus/cdna/Strongylocentrotus_purpuratus.Spur2.5.12.cdna.all.fa.gz
"tribolium_castaneum" => {'name' => 'Tribolium_castaneum','shortcut' => 'Tcas3', 'swcode' => 'TRICA', 'taxid' => '7070'},

"trichinella_spiralis" => {'name' => 'Trichinella_spiralis','shortcut' => 'Tspiralis1', 'swcode' => 'TRICA', 'taxid' => '7070'},
#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-15/fasta/trichinella_spiralis/pep/Trichinella_spiralis.Tspiralis1.15.pep.all.fa.gz
#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/tribolium_castaneum/cdna/Tribolium_castaneum.Tcas3.12.cdna.all.fa.gz
"trichoplax_adhaerens" => {'name' => 'Trichoplax_adhaerens','shortcut' => 'ASM15027v1', 'swcode' => 'TRIAD', 'taxid' => '10228'},
	#ftp://ftp.ensemblgenomes.org/pub/metazoa/release-12/fasta/trichoplax_adhaerens/cdna/Trichoplax_adhaerens.TRIAD1.12.cdna.all.fa.gz
	);
	print "\tdownloading metazoans\n";
	$current_ensemblDir = "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-$release/";
	foreach my $species_name(sort keys(%species2shortcut))
	{
		#'cdna' => 'ftp://ftp.ensembl.org/pub/release-66/fasta/ailuropoda_melanoleuca/cdna/Ailuropoda_melanoleuca.ailMel1.66.cdna.all.fa.gz',
		#'protein' => 'ftp://ftp.ensembl.org/pub/release-66/fasta/ailuropoda_melanoleuca/pep/Ailuropoda_melanoleuca.ailMel1.66.pep.all.fa.gz',
		#'gtf' => 'ftp://ftp.ensembl.org/pub/release-66/gtf/ailuropoda_melanoleuca/Ailuropoda_melanoleuca.ailMel1.66.gtf.gz',
		print "\tDownloading $species_name\n";
		#my $cndaSource = $species{$species_name}{'cdna'};
		my $speciesName2 = $species2shortcut{$species_name}{'name'};
		my $shortcut = $species2shortcut{$species_name}{'shortcut'};
		my $swcode = $species2shortcut{$species_name}{'swcode'};
		my $taxid = $species2shortcut{$species_name}{'taxid'};
		
		my $fileName = $speciesName2.".".$shortcut.".".$release.".".$cnda_file_ending;
		my $cndaSource = $current_ensemblDir."/fasta/".$species_name."/cdna/".$fileName; # remote
		my $cndaFile = $source_transcriptDir."/".$fileName; # local 

		unless(-e $cndaFile)
		{
			print "\tDownloading $cndaSource ($cndaFile)\n";
			getstore($cndaSource, $cndaFile);
		}	
		#my $proteinSource = $species{$species_name}{'protein'};
		$fileName = $speciesName2.".".$shortcut.".".$release.".".$protein_file_ending;
		my $proteinSource = $current_ensemblDir."/fasta/".$species_name."/pep/".$fileName; #remote
		my $proteinFile = $source_proteinDir."/".$fileName; # local 
		unless(-e $proteinFile)
		{
			print "\tDownloading $proteinSource\n";
			getstore($proteinSource, $proteinFile);
		}
		$fileName =~ s/\.gz//;
		$speciesName2toprint = $speciesName2;
		$speciesName2toprint =~ s/_/ /;
		$speciesListString .= "$speciesName2toprint\t$fileName\t$swcode\t$taxid\tENS\n";
		
		$fileName = $speciesName2.".".$shortcut.".".$release.".".$gtf_file_ending;
		#my $gtfFile = $source_gtfDir."/".basename($species{$species_name}{'gtf'});
		my $gtfFile = $source_gtfDir."/".$fileName;
		my $gtfSource = $current_ensemblDir."/gtf/".$species_name."/".$fileName;
		
		unless(-e $gtfFile){
			print "\tDownloading $gtfSource\n";
			getstore($gtfSource, $gtfFile);
		}
		
	}	
	write_to_file({file_name => $species_list, text => $speciesListString});
}
if($source eq 'jgi'){
	%species2shortcut = (
	"Monosiga_brevicollis" => {
		'proteins' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Monosiga_brevicollis/annotation/v1.0/Monbr1_all_proteins.fasta.gz',
		'transcripts' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Monosiga_brevicollis/annotation/v1.0/Monbr1_all_transcripts.fasta.gz',
		'gtf' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Monosiga_brevicollis/annotation/v1.0/Monbr1_all_models.gff.gz',
		'swcode' => 'MONBE',
		'taxid' => '81824',
	},
	"Branchiostoma_floridae" => {
		'proteins' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Branchiostoma_floridae/v1.0/proteins.Brafl1.fasta.gz',
		'transcripts' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Branchiostoma_floridae/v1.0/transcripts.Brafl1.fasta.gz',
		'gtf' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Bfloridae/v1.0/Brafl1.FilteredModels1.gff.gz',
		'swcode' => 'BRAFL',
		'taxid' => '7739',
	},
	"Capitella_sp" => {
		'proteins' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Capitella/v1.0/FilteredModelsv1.0.aa.fasta.gz',
		'transcripts' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Capitella/v1.0/FilteredModelsv1.0.na.fasta.gz',
		'gtf' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Capitella/v1.0/FilteredModelsv1.0.gff.gz',
		'swcode' => 'CAPIT',
		'taxid' => '73382',
	},
	"Helobdella_robusta" => {
		'proteins' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Helobdella_robusta/v1.0/proteins.Helro1_FilteredModels3.fasta.gz',
		'transcripts' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Helobdella_robusta/v1.0/transcripts.Helro1_FilteredModels3.fasta.gz',
		'gtf' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Helobdella_robusta/v1.0/Helobdella_robusta_FilteredModels3.gff.gz',
		'swcode' => 'HELRO',
		'taxid' => '6412'
	},
	
	);
	
	print "\tdownloading metazoans\n";
	$current_ensemblDir = "ftp://ftp.jgi-psf.org/pub/";
	foreach my $species_name(keys(%species2shortcut))
	{
		print "\tDownloading $species_name\n";
		#my $cndaSource = $species{$species_name}{'cdna'};
		my $speciesName2 = $species2shortcut{$species_name}{'name'};
		#my $shortcut = $species2shortcut{$species_name}{'shortcut'};
		my $swcode = $species2shortcut{$species_name}{'swcode'};
		my $taxid = $species2shortcut{$species_name}{'taxid'};
		my $fileName = basename($species2shortcut{$species_name}{'transcripts'});
		my $cndaSource = $species2shortcut{$species_name}{'transcripts'}; # remote
		my $cndaFile = $source_transcriptDir."/".$species_name.".cdna.all.fa.gz"; # local 
		unless(-e $cndaFile)
		{
			print "\tDownloading $cndaSource\n";
			getstore($cndaSource, $cndaFile);
			# gunzip file
			system("gunzip $cndaFile");
		}	
		#my $proteinSource = $species{$species_name}{'protein'};
		$fileName = basename($species2shortcut{$species_name}{'proteins'});

		my $proteinSource = $species2shortcut{$species_name}{'proteins'}; #remote
		my $proteinFile = $source_proteinDir."/".$species_name.".pep.all.fa.gz"; # local 
		unless(-e $proteinFile)
		{
			print "\tDownloading $proteinSource\n";
			getstore($proteinSource, $proteinFile);
			# gunzip file
			system("gunzip $proteinFile");
		}
		$fileName =~ s/\.gz//;
		$speciesName2toprint = $species_name;
		$speciesName2toprint =~ s/_/ /;
		
		$speciesListString .= "$speciesName2toprint\t$fileName\t$swcode\t$taxid\tJGI\n";
		
		$fileName = basename($species2shortcut{$species_name}{'gtf'});
		#my $gtfFile = $source_gtfDir."/".basename($species{$species_name}{'gtf'});
		my $gtfSource = $species2shortcut{$species_name}{'gtf'};
		my $gtfFile = $source_gtfDir."/".$species_name.".gtf.gz";
		
		unless(-e $gtfFile){
			print "\tDownloading $gtfSource\n";
			getstore($gtfSource, $gtfFile);
			# gunzip file
			system("gunzip $gtfFile");
		}
		
	}	
	write_to_file({file_name => $species_list, text => $speciesListString});
}
if($source eq 'ensembl'){
	
}



%species2shortcut = (
	"Monosiga_brevicollis" => {
		'proteins' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Monosiga_brevicollis/annotation/v1.0/Monbr1_all_proteins.fasta.gz',
		'transcripts' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Monosiga_brevicollis/annotation/v1.0/Monbr1_all_transcripts.fasta.gz',
		'gtf' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Monosiga_brevicollis/annotation/v1.0/Monbr1_all_models.gff.gz',
		'swcode' => 'MONBE',
		'taxid' => '81824',
	},
	"Branchiostoma_floridae" => {
		'proteins' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Branchiostoma_floridae/v1.0/proteins.Brafl1.fasta.gz',
		'transcripts' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Branchiostoma_floridae/v1.0/transcripts.Brafl1.fasta.gz',
		'gtf' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Bfloridae/v1.0/Brafl1.FilteredModels1.gff.gz',
		'swcode' => 'BRAFL',
		'taxid' => '7739',
	},
	"Capitella_sp" => {
		'proteins' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Capitella/v1.0/FilteredModelsv1.0.aa.fasta.gz',
		'transcripts' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Capitella/v1.0/FilteredModelsv1.0.na.fasta.gz',
		'gtf' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Capitella/v1.0/FilteredModelsv1.0.gff.gz',
		'swcode' => 'CAPIT',
		'taxid' => '73382',
	},
	"Helobdella_robusta" => {
		'proteins' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Helobdella_robusta/v1.0/proteins.Helro1_FilteredModels3.fasta.gz',
		'transcripts' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Helobdella_robusta/v1.0/transcripts.Helro1_FilteredModels3.fasta.gz',
		'gtf' => 'ftp://ftp.jgi-psf.org/pub/JGI_data/Helobdella_robusta/v1.0/Helobdella_robusta_FilteredModels3.gff.gz',
		'swcode' => 'HELRO',
		'taxid' => '6412'
	},
	
	);
# write table
# should contain columns: file name, ncbi tax id, species name
#print "printing species table\n";
my $species_file = "species_list.file";
my $print_string = '';
foreach my $species(keys(%species2shortcut)){
    my ($prot_file, $cds_file, $tax_id) = (basename($species2shortcut{$species}{'proteins'}),basename($species2shortcut{$species}{'transcripts'}),$species2shortcut{$species}{'taxid'}  );
    $prot_file =~ s/\.gz//;
    $cds_file =~ s/\.gz//;

                    $print_string .= "$warehouse_source_dir/proteins/$prot_file\t$warehouse_source_dir/transcripts/$cds_file\t$tax_id\t$species\n";
}
print $print_string;
#write_to_file({file_name => $species_file, text => $print_string});



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
