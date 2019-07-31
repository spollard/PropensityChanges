use strict;
use warnings;

use Data::Dumper;

use lib 'C:\Users\Stephen\Documents\Lab work\PerlModules';

use Sequences;


my $col_filename = "cols.txt";
my $fasta_filename = "hidden_states.txt";

open my $fi, "<", $col_filename;

my $names = [split " ", scalar(<$fi>)];
my $states = [split " ", scalar(<$fi>)];


my $seqs = {};
for my $i (0 .. $#$names){
    $seqs->{$names->[$i]} = $states->[$i]
}
print Dumper $seqs;

Sequences::ToFasta($seqs, $fasta_filename);


