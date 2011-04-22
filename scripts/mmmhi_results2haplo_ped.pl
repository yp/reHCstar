#!/usr/bin/perl -w


use strict;
use List::MoreUtils qw'mesh';
use List::MoreUtils qw'natatime';

use Getopt::Long;
use Pod::Usage;

my $help = 0;
my $pedigree_file= "";
my $bob_file= "";

GetOptions('pedigree-file|p=s' => \$pedigree_file,
           'bob-file|b=s' => \$bob_file,
           'help|?' => \$help
          ) or pod2usage(2);
pod2usage(1) if $help;


my %enc= ( '0 0' => '3',
           '1 1' => '0',
           '1 2' => '2',
           '2 2' => '1'
 );

if (${pedigree_file} =~ m/.gz$/) {
    open PED, "zcat ${pedigree_file} |" or
        die "Impossible to open pedigree file ${pedigree_file}!\n";
} else {
    open PED, "<${pedigree_file}" or
        die "Impossible to open pedigree file ${pedigree_file}!\n";
}
if (${bob_file} =~ m/.gz$/) {
    open BOB, "zcat ${bob_file} |" or
        die "Impossible to open haplotype configuration file ${bob_file}!\n";
} else {
    open BOB, "<${bob_file}" or
        die "Impossible to open haplotype configuration file ${bob_file}!\n";
}



my @indlist= ();
my %ped= ();
while (<PED>) {
    chomp;
    next if m/^\s*#/;
    last if m/^=====/;
    my @r= split(/\s/,$_,7);
    my @pedout= @r[0..5];
    my $id= $r[1];
    $ped{$id}= \@pedout;
    push(@indlist, $id);
}

<BOB>;

my %bobgens= ();
my $i= 1;
while (<BOB>) {
    chomp;
    last unless m/./;
    my @g= split(/  /,$_);
    $bobgens{"".$i}= \@g;
    ++$i;
}

my %phases= ();
while (<BOB>) {
    chomp;
    next unless m/^h_(\d+)_(\d+) = ([0-1])$/;
    $phases{(0+$1)."-".(0+$2)}= 0+$3;
}

my %haplotypes= ();
foreach my $id (sort {$a <=> $b} keys %ped) {
    my $f= $ped{$id}->[2];
    my $m= $ped{$id}->[3];
    next unless ($f == 0) && ($m == 0);
    my @hp= ();
    my @hm= ();
    foreach my $g (@{$bobgens{$id}}) {
        if ($g =~ m/([01])([01])/) {
            push(@hp, 1+$1);
            push(@hm, 1+$2);
        } elsif ($g =~ m/([01])\?/) {
            push(@hp, 1+$1);
            push(@hm, 1);
        } elsif ($g =~ m/\?([01])/) {
            push(@hp, 1);
            push(@hm, 1+$1);
        } elsif ($g =~ m/\?\?/) {
            push(@hp, 1);
            push(@hm, 1);
        } elsif ($g =~ m/\*\*/) {
            push(@hp, 1);
            push(@hm, 2);
        } else {
            die "Invalid genotype $g\n";
        }
    }
    $haplotypes{$id}= [ \@hp, \@hm ];
}
while ((scalar keys %haplotypes) < (scalar keys %ped)) {
    foreach my $id (sort {$a <=> $b} keys %ped) {
        my $f= $ped{$id}->[2];
        my $m= $ped{$id}->[3];
        next unless defined $haplotypes{$f} && defined $haplotypes{$m};
        my @hp= ();
        my @hm= ();
        my $l= 0;
        foreach my $g (@{$bobgens{$id}}) {
            if ($g =~ m/([01])([01])/) {
                push(@hp, 1+$1);
                push(@hm, 1+$2);
            } elsif ($g =~ m/([01])\?/) {
                push(@hp, 1+$1);
                push(@hm, $haplotypes{$m}->[$phases{"$m-$id"}]->[$l]);
            } elsif ($g =~ m/\?([01])/) {
                push(@hp, $haplotypes{$f}->[$phases{"$f-$id"}]->[$l]);
                push(@hm, 1+$1);
            } elsif ($g =~ m/[\?\*][\?\*]/) {
                push(@hp, $haplotypes{$f}->[$phases{"$f-$id"}]->[$l]);
                push(@hm, $haplotypes{$m}->[$phases{"$m-$id"}]->[$l]);
            } else {
                die "Invalid genotype $g\n";
            }
            ++$l;
        }
        $haplotypes{$id}= [ \@hp, \@hm ];
    }
}

foreach my $id (@indlist) {
    print(join("\t", @{$ped{$id}})."\t");
    my $hr= $haplotypes{$id};
    my $it= natatime(2, mesh(@{$hr->[0]}, @{$hr->[1]}));
    my $i= 0;
    while (my @h= $it->()) {
        print("\t") if ($i>0);
        print(join("|", @h));
        ++$i;
    }
    print("\n");
}


__END__

=head1 SYNOPSIS

mmmhi_results2haplo_ped.pl [--pedigree-file|-p]=FILE [--bob-file|-b]=FILE

 Read from 'pedigree-file' the PED-format pedigree structure and
 combine it with the haplotype configuration read from the file 'bob-file'
 in the Bob format.

=cut
