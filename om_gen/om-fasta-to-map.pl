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
        "   -c, --chroms <list>        Chromosomes to include, \"-\" for all\n",
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
        chroms => [ qw(1 10 11 12 13 14 15 16 17 18 19 2 3 4 5 6 7 8 9 MT X Y) ],
    };
    for (my $i=0; $i<@{$$opts{chroms}}; $i++) { $$opts{chroms}[$i] = "chr".$$opts{chroms}[$i]; }
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( $arg eq '-c' || $arg eq '--chroms' ) { $$opts{chroms} = [split(/,/,shift(@ARGV))]; next; }
        if ( $arg eq '-f' || $arg eq '--fasta' ) { $$opts{fasta} = shift(@ARGV); next; }
        if ( $arg eq '-m' || $arg eq '--motif' ) { $$opts{motif} = shift(@ARGV); next; }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( ! -e $$opts{fasta} ) { error("Uh: no such file: $$opts{fasta}\n"); }
    return $opts;
}

sub rev_strand
{
    my ($seq) = @_;
    my %map = (a=>'t',c=>'g',g=>'c',t=>'a',u=>'a');
    my $out;
    my $len = length($seq);
    for (my $i=$len-1; $i>=0; $i--)
    {
        my $base = substr($seq,$i,1);
        $out .= exists($map{$base}) ? $map{$base} : $base;
    }
    return $out;
}

sub read_fasta
{
    my ($opts) = @_;
    my $all_chrs = 0;
    my %chrs = map { $_ => 1 } @{$$opts{chroms}};
    if ( exists($$opts{chroms}) && $$opts{chroms}[0] eq '-' ) { $all_chrs = 1; }
    my $fwd_motif = lc($$opts{motif});
    my $rev_motif = rev_strand($fwd_motif); #"ccatgg"
    my $motif_len = length($fwd_motif);
    my $skip  = 0;
    my $off   = 0;
    my $buf   = '';
    my $chr   = '';
    my $prev_hit = 0;
    my $cmd = $$opts{fasta}=~/gz$/ ? "gunzip -c $$opts{fasta} |" : "<$$opts{fasta}";
    open(my $fh,$cmd) or error("$cmd: $!");
    while (my $line=<$fh>)
    {
        $line = lc($line);
        if ( $line=~/^>(\S+)/ )
        {
            if ( $all_chrs or exists($chrs{$1}) ) { $skip = 0; }
            else { $skip = 1; }
            $off = 0;
            $buf = '';
            $prev_hit = 0;
            $chr = $1;
            next;
        }
        if ( $skip ) { next; }
        chomp($line);
        $buf .= $line;
        my $fwd_idx = -1;
        my $rev_idx = -1;
        my $idx     = 0;
        while (1)
        {
            $fwd_idx = index($buf,$fwd_motif,$idx);
            $rev_idx = index($buf,$rev_motif,$idx);

            if ( $fwd_idx==-1 && $rev_idx==-1 ) { last; }   # no index

            if ( $fwd_idx==-1 ) { $fwd_idx = 1e10; }
            if ( $rev_idx==-1 ) { $rev_idx = 1e10; }
            $idx = $fwd_idx < $rev_idx ? $fwd_idx : $rev_idx;

            my $hit = $off + $idx;
            $idx += $motif_len;
            if ( $prev_hit > 0 ) # don't output the first fragment, the start of the fasta reference is arbitrary
            {
                #if ($hit-$prev_hit > 2000)  {
                printf "$chr\t$prev_hit\t$hit\t%.3f\n", ($hit-$prev_hit-$motif_len)*1e-3;
                #}
            }
            $prev_hit = $hit;

        }

        my $len = length($buf) - $motif_len + 1;
        substr($buf,0,$len,'');
        $off += $len;
        $idx  = 0;
    }
    close($fh);
}



