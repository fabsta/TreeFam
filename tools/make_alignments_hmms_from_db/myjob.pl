#!/usr/bin/perl 
#===============================================================================
#
#         FILE: DumpPreviousTreeFamAlnHMMs.pl
#
#        USAGE: ./DumpPreviousTreeFamAlnHMMs.pl
#
#  DESCRIPTION: Dumps all alns and hmms from previous release
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Fabian Schreiber (), fs9@sanger.ac.uk
# ORGANIZATION:
#      VERSION: 1.0
#      CREATED: 06/07/2012 01:13:18 PM
#     REVISION: ---
#===============================================================================

use strict;
use File::Temp qw/ tempfile tempdir /;
use warnings;
use TreeFam::DBConnection;
use TreeFam::FamilyHandle;
use Bio::SeqIO;
use File::Basename;
use Data::Dumper;
# Iterate over all families

my @all_jobs = `cat treefam8_allfamilies3.txt`; 
my $treefamFamilyID_number       = $ARGV[0];
$treefamFamilyID_number =~ s/input\.//;
#print "will try $treefamFamilyID_number\n";
my $treefamFamilyID = $all_jobs[$treefamFamilyID_number-1];
chomp($treefamFamilyID);
#print "will try $treefamFamilyID\n";

my $worker_temp_directory = "/tmp/treefam/$treefamFamilyID";
`mkdir -p $worker_temp_directory`;
#print "looking at family $treefamFamilyID\n";

#($fh, $alignmentFile) = tempfile();
my $tmpalignmentFile     = $worker_temp_directory . "/" . $treefamFamilyID . ".tmpaln";
my $alignmentFile        = $worker_temp_directory . "/" . $treefamFamilyID . ".aln";
my $stockholmFile        = $worker_temp_directory . "/" . $treefamFamilyID . ".sto";
#my $hmmFile              = $worker_temp_directory . "/" . $treefamFamilyID . ".hmm";
my $hmmFile              = $worker_temp_directory . "/hmmer.hmm";
my $realignment_required = -1;

#iterate over all treefam version -1 things
#print STDERR "Getting clean alignment for $treefamFamilyID...\n";
my $dbc  = new TreeFam::DBConnection();
my $trh = $dbc->get_TreeHandle();
#my $famh = $dbc->get_FamilyHandle();
my $fam_type = $treefamFamilyID =~ /^TF1/ ? "SEED":"CLEAN";

if ( defined( my $tree = $trh->get_by_id($treefamFamilyID, $fam_type) ) )
{
    if(-e $hmmFile && -s $hmmFile){
        die "Family already exists. Skipping\n";

        }
    #my @ListOfSequences = $famh->get_gene_sequences($treefamFamilyID);
    my @ListOfLeaves = $tree->get_leaves;
    my @ListOfSequences;
    foreach my $leaf (@ListOfLeaves){
        #print "leaf: ".$leaf->name." id: ".$leaf->id." symbol: ".$leaf->sequence_id." sequence: ".$leaf->sequence."\n";
        next if !$leaf || $leaf eq '';
        next if !(eval { $leaf->can("sequence") });
        my $seq = $leaf->sequence;
        next if !$seq || $seq eq '';
        my $id = $leaf->sequence_id;
        push(@ListOfSequences, Bio::Seq->new( -seq => $seq,-id  => $id));    
    }
    if ( !@ListOfSequences )
    {
        print STDERR "\tCould not get gene sequences for family $treefamFamilyID\n";
    }
    my $seqout = Bio::SeqIO->new( -format => 'Fasta', -file => ">$alignmentFile" );
    my $length;
    foreach (@ListOfSequences)
    {
        $length = length( $_->seq ) if !defined($length);
        $seqout->write_seq($_);

        if ( $length != length( $_->seq ) )
        {
            #print "Sequences not all same length\n";
            $realignment_required = 1;
        }
    }
    if ( $realignment_required == 1 )
    {
        #print "\tseq length mismatch, error in cigar string. Realignment required.\n";

        #TODO replace with module
        my $muscle_cmd = "muscle -quiet -in $alignmentFile -out $tmpalignmentFile ";
        system($muscle_cmd);
        if ( !-e $tmpalignmentFile || !-s $tmpalignmentFile )
        {
            die "Could not compute alignment using $muscle_cmd \n";
        }
        `mv $tmpalignmentFile $alignmentFile`;
    
    }
    print "STATISTICS:$treefamFamilyID\t".scalar(@ListOfSequences)."\n";
    &fasta2stockholm( $alignmentFile, $stockholmFile );
    if ( !-e $stockholmFile || !-s $stockholmFile )
    {
        die "\tcould not convert to stockholm file $stockholmFile\n";
    }
    ## do hmmbuild
    print "building HMM\n";
   
    `rm $hmmFile` if -e $hmmFile;
    my $hmmbuild = "/software/ensembl/compara/hmmer-2.3.2/src/hmmbuild";
    #my $hmmbuild = "/Users/fs9/bin/source/hmmer-2.3.2/src/hmmbuild";
    my $cmd      = $hmmbuild;
    $cmd .= ' --amino ';
    $cmd .= $hmmFile;

    $cmd .= " " . $stockholmFile;
    $cmd .= " 2>&1 > /dev/null";

    #print("$cmd\n");
    system($cmd);
    my $lustre_dir = "/lustre/scratch109/sanger/fs9/treefam8_hmms_dir/books";
    #print "copying files to lustre: $lustre_dir\n";
    my $cp_cmd = "cp -r ".dirname($hmmFile)." $lustre_dir";
    #print ("$cp_cmd\n");
    system($cp_cmd);
}
else
{
    print STDERR "WARNING: family $treefamFamilyID is not in the database\n";
}

sub fasta2stockholm
{
    my $fasta_file    = shift;
    my $stockholmFile = shift;

    my %seq;
    my @name;
    my $name;
    open FASTA, "<$fasta_file" or die "Couldn't open '$fasta_file': $!";
    while (<FASTA>)
    {
        if (/^\s*>\s*(\S+)/)
        {
            $name = $1;
            die "Duplicate name: $name" if defined $seq{$name};
            push @name, $name;
        }
        else
        {
            if ( /\S/ && !defined $name )
            {
                warn "Ignoring: $_";
            }
            else
            {
                s/\s//g;
                $seq{$name} .= $_;
            }
        }
    }
    close FASTA;

    # print Stockholm output
    open my $STOCKHOLM_FILE, ">$stockholmFile"
        or die "Couldn't open '$stockholmFile': $!";
    print {$STOCKHOLM_FILE} "# STOCKHOLM 1.0\n";
    foreach my $name (@name)
    {
        print {$STOCKHOLM_FILE} $name, " ", $seq{$name}, "\n";
    }
    print {$STOCKHOLM_FILE} "//\n";

    close $STOCKHOLM_FILE;
    return ( -e $stockholmFile && -s $stockholmFile ) ? $stockholmFile : undef;

}


#sub get_leaf_sequences{
    #my ($tree) = (@_);
    #print "looking at tree $tree\n";

    #use Bio::Phylo::IO;
     #my $string = '((A,B),C);';
     #my $forest = Bio::Phylo::IO->parse(
        #-format => 'newick',
        #-string => $string
     #);
     #my $tree = $forest->first;
    #my @terminals = @{ $tree->get_terminals };
#}
