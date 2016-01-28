#!/usr/bin/env perl
#
# Author: petr.danecek@sanger
#

use strict;
use warnings;
use Carp;

my $opts = parse_params();
read_fasta($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print 
        "Usage: om-fasta-to-map [OPTIONS]\n",
        "Options:\n",
        "   -f, --fasta <file.fa>      \n",
        "   -m, --motif <string>       \n",
        "   -h, -?, --help             This help message.\n",
        "\n";
    exit -1;
}


sub parse_params
{
    my $opts = 
    {
        fasta => '/lustre/scratch107/vrpipe/refs/mouse/ncbim37/NCBIM37_um.fa',
        motif => 'GGTACC',
        chroms => [ qw(chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrMT chrX chrY 1 10 11 12 13 14 15 16 17 18 19 2 3 4 5 6 7 8 9 MT X Y) ],
    };
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( $arg eq '-f' || $arg eq '--fasta' ) { $$opts{fasta} = shift(@ARGV); next; }
        if ( $arg eq '-m' || $arg eq '--motif' ) { $$opts{motif} = shift(@ARGV); next; }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( ! -e $$opts{fasta} ) { error("Uh: no argument for -f?\n"); }
    return $opts;
}

sub read_fasta
{
    my ($opts) = @_;
    my %chrs = map { $_ => 1 } @{$$opts{chroms}};
    my $motif = $$opts{motif};
    my $skip  = 0;
    my $off   = 0;
    my $buf   = '';
    my $chr   = '';
    my $prev_hit = undef;
    open(my $fh,'<',$$opts{fasta}) or error("$$opts{fasta}: $!");
    while (my $line=<$fh>)
    {
        if ( $line=~/^>(\S+)/ )
        {
            if ( exists($chrs{$1}) ) { $skip = 0; }
            else { $skip = 1; }
            $off = 0;
            $buf = '';
            $prev_hit = undef;
            $chr = $1;
        }
        if ( $skip ) { next; }
        chomp($line);
        $buf .= $line;
        my $idx = -1;
        do
        {
        	my $search_offset = 0; #using only "$idx+length($motif)" misses motifs at the beginnings of lines
        	if ($idx >=0 ){ $search_offset = $idx+length($motif);}  
            $idx = index($buf,$motif,$search_offset);
            if ( $idx!=-1 )
            {
                my $hit = $off + $idx;
                if ( defined $prev_hit )
                {
                    printf "$chr\t$prev_hit\t$hit\t%.3f\n", ($hit-$prev_hit-length($motif))*1e-3;
                }
                $prev_hit = $hit;
            }
        }
        while ($idx!=-1);

        my $len = length($line);
        substr($buf,0,$len,'');
        $off += $len;
    }
    close($fh);
}



