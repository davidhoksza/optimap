#!/usr/bin/env perl
#
# Author: petr.danecek@sanger
#

use strict;
use warnings;
use Carp;
#use Random;

my $opts = parse_params();
sim_data($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print 
        "Usage: om-sim-data [OPTIONS]\n",
        "Options:\n",
        "   -m, --map <ref.map.gz>     Reference map, output of om-fasta-to-map\n",
        "   -d, --digestion n          Digestion rate (0;1]\n",
        "   -h, -?, --help             This help message.\n",
        "\n";
    exit -1;
}


sub parse_params
{
    my $opts = 
    {
        dgst_rate  => 0.8,     # fraction of sites digested
        brk_rate   => 1e-5,    # random breaks for the Poisson process, rare event
        min_nsites => 12,      # minimum number of fragments per molecule
        min_mlen   => 50,     # minimum molecule length [kb]
        max_mlen   => 250,     # maximum molecule length [kb]
        map_fname  => 'ref.map.gz', # reference map
    };
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( $arg eq '-m' || $arg eq '--map' ) { $$opts{map_fname} = shift(@ARGV); next; }
        if ( $arg eq '-d' || $arg eq '--digestion' ) { $$opts{dgst_rate} = shift(@ARGV); next; }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( ! -e $$opts{map_fname} ) { error("Uh: $$opts{map_fname}\n"); }
    return $opts;
}

sub next_break
{
    my ($rate) = @_;
    return int(-log(1.0 - rand())/$rate);
}

sub sim_data
{
	#srand(1);
    my ($opts) = @_;
    my $brk_at = next_break($$opts{brk_rate});
    open(my $fh,"gunzip -c $$opts{map_fname} |") or error("$$opts{map_fname}: $!");
    while (1)
    {
    	#my $brk_at = next_break($$opts{brk_rate});
    	
        my $mlen = $$opts{min_mlen} + int(rand($$opts{max_mlen} - $$opts{min_mlen}));   # random molecule length
        my @buf = ();
        my @pos = ();
        my @pos_from = ();
        my $len = 0;
        my $prev_chr = undef;
        my $line;
        while ($line=<$fh>)
        {
            chomp($line);
            my ($chr,$from,$to,$ilen) = split(/\t/,$line);
            if ( defined $prev_chr && $prev_chr ne $chr ) { last; }
            push @buf,$ilen;
            push @pos,$to;
            push @pos_from,$from;
            $len += $ilen;
            if ( $len >= $mlen ) { last; }
            $prev_chr = $chr;
        }
        if ( !defined $line ) { last; }
        if ( @buf<=$$opts{min_nsites} ) { next; }       # in real data we observed at least min_nsites per molecule        
        my $pos = "debug\t$prev_chr\t$pos_from[0]\t";
        my $out = "\tKpnI\tK\t";
        $len = 0;                   # uncleaved length [kbp]
        my $prn_len = 0;            # molecule length printed so far [bp]
        my $missed_sites = 0;
        for (my $i=0; $i<@buf; $i++)
        {
            my $xkb = $buf[$i];       # length of the current reference fragment
            my $xb  = 1000*$xkb;
            if ( $$opts{brk_rate} && $prn_len+$xb+$len*1000 > $brk_at )    # time for next random break
            {
                $pos .= " b:".($pos_from[0]+$brk_at);     # approximate absolute position of the random break; for debugging
                $out .= "\t".(($brk_at - $prn_len)/1000);               
                
                           
            	$prn_len += $xb + $len*1000;
                
                $pos .= " e:$pos[$i]";      # absolute position of the enzymatic cleavage; for debugging
            	$out .= "\t".(($prn_len - $brk_at)/1000);
            	$len = 0;
            	
            	while ($brk_at < $prn_len) {$brk_at += next_break($$opts{brk_rate});}
            	#$brk_at = next_break($$opts{brk_rate});
                next;
            }
            if ( rand(1.0)>$$opts{dgst_rate} && $missed_sites < 4 ) # random cleavage site omission
            { 
                $missed_sites += 1;
                $len += $xkb;
                next; 
            }
            $pos .= " e:$pos[$i]";      # absolute position of the enzymatic cleavage; for debugging
            $out .= "\t".($xkb+$len);
            
            $prn_len += $xb+$len*1000;
            $len  = 0;
        }
        $brk_at -= $len;
        print $pos,"\n";        # for debugging, not part of the real data
        print $out,"\n";
    }
    close($fh);
}



