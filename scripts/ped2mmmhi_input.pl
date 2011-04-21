#!/usr/bin/perl -w


use strict;

my %genders= ( '1' => 'M', '2' => 'F' );
my %enc= ( '0 0' => '3',
           '1 1' => '0',
           '1 2' => '2',
           '2 2' => '1'
 );

my %ped= ();
my %gen= ();
my $glen= 0;
while (<>) {
    chomp;
    next if m/^\s*#/;
    my @r= split(/\s/,$_,7);
    my @pedout;
    my $id= $r[1];
    push(@pedout, $id);
    push(@pedout, $genders{$r[4]});
    if ($r[2] ne '0') {
        push(@pedout, $r[2]);
    } else {
        push(@pedout, 'N');
    }
    if ($r[3] ne 0) {
        push(@pedout, $r[3]);
    } else {
        push(@pedout, 'N');
    }
    $ped{$id}= \@pedout;

    my @g= split(/\t/,$r[6]);
    my $newg= "";
    foreach my $s (@g) {
        $newg .= $enc{$s}." ";
    }
    chomp($newg);
    $gen{$id}= $newg;
    if ($glen == 0) {
        my $g= $newg;
        $g =~ s/ //g;
        $glen= length($g);
    }
}

print(scalar(keys %ped)."\n");
foreach my $id (sort {$a <=> $b} keys %ped) {
    print(join("\t",@{$ped{$id}})."\n");
}

print(STDERR scalar(keys %gen)." ".$glen."\n");
foreach my $id (sort {$a <=> $b} keys %gen) {
    print(STDERR $gen{$id}."\n");
}
