#!/usr/bin/perl

use strict;
use warnings;

use Bio::Phylo::IO;
use JSON;

my $forest=Bio::Phylo::IO->parse(-file=>"example.newick", -format=>"newick");

while (my $tree=$forest->next) {
    my $out=[];
    my $children=$out;
    my $cur;
    my $parent;
    $tree->visit_breadth_first(
                               -pre=>sub { $cur={ name=>shift->get_name }; push @$children, $cur; },
                               -post =>sub { $cur->{children}=[]; $parent=$cur;  $children=$cur->{children} },
                              );
    print JSON->new->pretty->encode($out);
}

