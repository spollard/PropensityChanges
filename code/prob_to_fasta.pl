use strict;
use warnings;

use Data::Dumper;

use lib 'C:\Users\Stephen\Documents\Lab work\PerlModules';

use Sequences;


my $probs_filename = $ARGV[0] // "probs.txt";
my $fasta_filename = $ARGV[1] // "rounded_hidden_states.txt";


my $seqs = Sequences::FromFasta($probs_filename);

for my $name (keys(%$seqs)) {
    $seqs->{$name} = $seqs->{$name} >= 0.5 ? 1 : 0;
    print "$name\t$seqs->{$name}\n";
}

Sequences::ToFasta($seqs, $fasta_filename);


