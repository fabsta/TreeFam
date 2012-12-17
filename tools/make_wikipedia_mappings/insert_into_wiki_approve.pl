#!/usr/bin/perl
use Data::Dumper;
use Scalar::Util;
use Getopt::Long;
use strict;
use warnings;
use JSON;
use File::Temp qw/ tempdir tempfile /;
use Log::Log4perl qw(get_logger :levels);
use Config::General;

BEGIN {
      my $pfad = "/nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/PfamSchemata/";
        push(@INC,$pfad);
}

use WikiApprove;
my $config_file = 'wiki.conf';
# set up logging
my $logger_conf = q(
  log4perl.logger                   = INFO, Screen
  log4perl.appender.Screen          = Log::Log4perl::Appender::Screen
  log4perl.appender.Screen.layout   = Log::Log4perl::Layout::PatternLayout
  log4perl.appender.Screen.layout.ConversionPattern = %M:%L %p: %m%n
);

Log::Log4perl->init( \$logger_conf );

my $log = get_logger();
# parse the config and get the section relevant to the web_user database
my $cg        = Config::General->new($config_file);
my %config    = $cg->getall;
my $wa_conf   = $config{wiki_approve};

#print Dumper %config;
print Dumper $wa_conf;
#exit;

# get all of the database connections that we'll need
my $wa_schema = 
  WikiApprove->connect( 
    "dbi:mysql:$wa_conf->{db_name}:$wa_conf->{db_host}:$wa_conf->{db_port}", 
    $wa_conf->{username},
    $wa_conf->{password}
  );



my $debug = 0;
my $clean_old_entries  = 1;

my $mappings_dir = "/nfs/users/nfs_f/fs9/bin/perl_modules/ensembl_main/treefam_tools/make_wikipedia_mappings/";
my $mappings_file = "$mappings_dir/mappings.txt";

print "Check existence of mapping files" if $debug;
die "File desc missing\n" if(!-e $mappings_file|| !-s $mappings_file);
my %treefam_map = read_mappings($mappings_file);
die "problem reading $mappings_file file\n" if !keys(%treefam_map);
print "well, we found ".keys(%treefam_map)." mappings\n";
#print Dumper %treefam_map;

#exit;
my $counter = 0;
foreach my $acc ( keys %treefam_map ) { 
  my $titles = $treefam_map{$acc};

  foreach my $title ( @$titles ) {
    $log->debug( "checking TreeFam entry/title: |$acc|$title|" );
    eval {
        print "add_row( $acc, $title, 'treefam' )\n";
        add_row( $acc, $title, 'treefam' );
    };
    if ( $@ ) {
      $log->logwarn( $@ );
    }
    #last if $counter++ > 2;
  }
}

#inserting into gene_tree_root_tag

#-------------------------------------------------------------------------------
#- functions -------------------------------------------------------------------
#-------------------------------------------------------------------------------

# adds a row to the "article_mapping" table with the specified accession, title
# and "db" values, and a corresponding row to "wikipedia". All rows in the 
# "article_mapping" table with the given accession are deleted before an attempt
# to add the new ones. If the delete fails, a warning is issues. If inserting 
# into either "article_mapping" or "wikipedia" fails, an exception is thrown.

sub add_row {
  my ( $acc, $title, $db ) = @_;

  # this should be in a transaction

  $wa_schema->resultset('ArticleMapping')
            ->update_or_create( { title     => $title,
                                  accession => $acc,
                                  db        => $db },
                                { key => 'primary' } )
    or die "error: failed to add mapping for '$acc' --> '$title' -> $db";
                                  
# TODO updating a treefam article will set both pfam_status and rfam_status = "inactive"
    print "update $title set treefam_status = 'active'\n"; 
  $wa_schema->resultset('Wikipedia')
            ->update_or_create( { title       => $title,
                                  approved_by => 'new',
                                  #pfam_status => $db eq 'pfam' ? 'active' : 'inactive', 
                                  #rfam_status => $db eq 'rfam' ? 'active' : 'inactive' },
                                  treefam_status => 'active'},
                                { key => 'primary' } )
  or die "error: failed to add wikipedia row for '$acc', '$title'";
}

sub read_mappings{
	my ($file) = (@_);
	my %entrez2name;
    my @lines_of_file = `cat $file`;
    my $counter = 0;
	foreach my $line(@lines_of_file){
		chomp($line);
		my ($tf,$wiki) = split(/\t/,$line);
	    if(!defined($tf) || $tf eq '' || !defined($wiki) || $wiki eq ""){
            warn "empty something reading tf:$tf and wiki:$wiki\n";
        }
        else{
            push(@{$entrez2name{$tf}},$wiki);
            #last if $counter++ > 3;
        }
}
	#print "Found ".keys(%entrez2name)." entrez infos\n";
    return %entrez2name;
}
