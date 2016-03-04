#!/usr/bin/env perl
#
# Author: petr.danecek@sanger
#
#   bgzip *.fasta *bed
#   samtools faidx C57B6J_1504.fasta.gz
#   samtools faidx CASTEiJ_1504.fasta.gz
#   zcat CASTEiJ_1504-C57B6J_chr2.bed.gz | sort -k1,1d -k2,2n -k4,4n | bgzip -c > CASTEiJ_1504-C57B6J_chr2.dst-sorted.bed.gz
#

use strict;
use warnings;
use Carp;
use FaSlice;

my $opts = parse_params();
process($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print 
        "Usage: script [OPTIONS]\n",
        "Options:\n",
        "   -b, --bed <file>         Lifted regions.\n",
        "   -c, --chr <string>       Current chromosome.\n",
        "   -f, --fasta <file>       Source fasta file.\n",
        "   -h, -?, --help           This help message.\n",
        "\n";
    exit -1;
}
sub parse_params
{
    my $opts = 
    {
        bed       => 'CASTEiJ_1504-C57B6J_chr2.dst-sorted.bed.gz',
        src_fasta => 'CASTEiJ_1504.fasta.gz',
        dst_fasta => 'C57B6J_1504.fasta.gz',
        chr       => '2',
        motif     => 'GGTACC',
        max_gap   => 500,     # as long as it does not contain the restriction site [bp]
        min_len   => 1e6,     # min block to report [bp]
    };
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( $arg eq '-c' || $arg eq '--chr' ) { $$opts{chr}=shift(@ARGV); next }
        if ( $arg eq '-b' || $arg eq '--bed' ) { $$opts{bed}=shift(@ARGV); next }
        if ( $arg eq '-f' || $arg eq '--fasta' ) { $$opts{fasta}=shift(@ARGV); next }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    $$opts{src_fa} = FaSlice->new(file=>$$opts{src_fasta});
    $$opts{dst_fa} = FaSlice->new(file=>$$opts{dst_fasta});
    return $opts;
}

sub process
{
    my ($opts) = @_;
    print "# [1]CAST_beg\t[2]CAST_end\t[3]C57BL_beg\t[4]C57BL_end\t[5]CAST_len (kb)\t[6]C57BL_len (kb)\n";
    open(my $bed,"zcat $$opts{bed} |") or error("zcat $$opts{bed}: $!");
    while (my $line=<$bed>)
    {
        my ($chr,$dst_beg,$dst_end,$src) = split(/\s+/,$line);
        chomp($src);
        my ($src_beg,$src_end) = split(/-/,$src);
        if ( $chr ne $$opts{chr} ) { next; }
        flush_buffer($opts,{dst_beg=>$dst_beg,dst_end=>$dst_end,src_beg=>$src_beg,src_end=>$src_end});
    }
    close($bed) or error("close failed: zcat $$opts{bed}");

    for my $key (qw(tot outlier adjacent skip1 skip2 skip3 good_enough out))
    {
        printf STDERR "%.2f \t ${key}_src\n", $$opts{"${key}_src"}*1e-6;
        printf STDERR "%.2f \t ${key}_dst\n", $$opts{"${key}_dst"}*1e-6;
    }
}

sub has_motif
{
    my ($opts,$fa,$beg,$end) = @_;
    if ( $beg >= $end ) { return 0; }
    $beg -= length($$opts{motif});
    $end += length($$opts{motif});
    if ( $beg<0 ) { $beg = 0; }
    my $seq = $fa->get_slice($$opts{chr},$beg,$end);
    if ( index($seq,$$opts{motif})>=0 ) { return 1; }
    return 0;
}

# Relying on the records being sorted by $$reg{dst_beg},$$reg{src_beg}
sub flush_buffer
{
    my ($opts,$reg) = @_;

    $$opts{tot_src} += $$reg{src_end} - $$reg{src_beg};
    $$opts{tot_dst} += $$reg{dst_end} - $$reg{dst_beg};

    if ( !defined $$opts{prev} ) { $$opts{prev} = { %$reg }; return; }
    my $prev = $$opts{prev};

    if ( ($$reg{src_beg} - $$prev{src_beg}) > 1e6 )
    { 
        # outlier, skip
        $$opts{outlier_src} += $$reg{src_end} - $$reg{src_beg};
        $$opts{outlier_dst} += $$reg{dst_end} - $$reg{dst_beg};
        return; 
    }
    
    if ( $$reg{src_beg}==$$prev{src_end} )
    {
        if ( $$reg{dst_beg}==$$prev{dst_end} or !has_motif($opts,$$opts{dst_fa},$$prev{dst_end},$$reg{dst_beg}) )
        {
            # well behaved region: perfect translation or a deletion/translocation w/o restriction
            # site in cast
            $$prev{src_end} = $$reg{src_end}; 
            $$prev{dst_end} = $$reg{dst_end}; 
            $$opts{adjacent_src} += $$reg{src_end} - $$reg{src_beg};
            $$opts{adjacent_dst} += $$reg{dst_end} - $$reg{dst_beg};
        }
        else
        {
            $$opts{skip1_src} += $$reg{src_end} - $$reg{src_beg};
            $$opts{skip1_dst} += $$reg{dst_end} - $$reg{dst_beg};
        }
        return;
    }
    if ( $$prev{src_end} > $$reg{src_beg} )
    {
        # misplaced, skip
        $$opts{skip2_src} += $$reg{src_end} - $$reg{src_beg};
        $$opts{skip2_dst} += $$reg{dst_end} - $$reg{dst_beg};
        return;
    }

    if ( $$reg{src_beg} - $$prev{src_end} < $$opts{max_gap} && $$reg{dst_beg} - $$prev{dst_end} < $$opts{max_gap} )
    {
        my $good_enough = 0;

        # small gap
        if ( !has_motif($opts,$$opts{src_fa},$$prev{src_end},$$reg{src_beg}) ) { $good_enough++; }
        if ( !has_motif($opts,$$opts{dst_fa},$$prev{dst_end},$$reg{dst_beg}) ) { $good_enough++; }
        if ( $good_enough == 2 )
        {
            $$prev{src_end} = $$reg{src_end}; 
            $$prev{dst_end} = $$reg{dst_end}; 
            $$opts{good_enough_src} += $$reg{src_end} - $$reg{src_beg};
            $$opts{good_enough_dst} += $$reg{dst_end} - $$reg{dst_beg};
        }
        else
        {
            $$opts{skip3_src} += $$reg{src_end} - $$reg{src_beg};
            $$opts{skip3_dst} += $$reg{dst_end} - $$reg{dst_beg};
        }
        return;
    }

    printf "%d\t%d\t%d\t%d\t%.1f\t%.1f\n", 
        $$prev{src_beg},$$prev{src_end},$$prev{dst_beg},$$prev{dst_end},
        1e-3*($$prev{src_end}-$$prev{src_beg}),1e-3*($$prev{dst_end}-$$prev{dst_beg});

    $$opts{out_src} += $$prev{src_end} - $$prev{src_beg};
    $$opts{out_dst} += $$prev{dst_end} - $$prev{dst_beg};

    $$opts{prev} = { %$reg };
}



