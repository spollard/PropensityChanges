our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;



use Sequences;

my $in = "../data/2seq";
my $out = "../data/100seq";

my $sites = [1 .. 100];

my $seqs = Sequences::FromFasta($in);
my $seqs = Sequences::Select($seqs, $sites);
Sequences::ToFasta($seqs, $out);
