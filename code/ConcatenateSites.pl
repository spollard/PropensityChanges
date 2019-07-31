use strict;
use warnings;


our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;

use Sequences;

unless (caller) {
    my $sites_dir = $ARGV[0] // "../data/simulated_100_sites_100_taxa_4/sites/" ;

    my $sequences_filename = $ARGV[1] // $sites_dir . "sequences.fasta";
    
    #~ my $fastas = [<$sites_dir*sequences.fasta>];
    #~ print join "\n", @$fastas;

    opendir(my $dir, $sites_dir) or die $!;
    my $site_dirs = [grep { -d $_} map {"$sites_dir/$_"} grep {/^\d+/} readdir($dir)];
    my $site_files = [map {$_ . "/sequences"} sort {($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0]} @$site_dirs];

    
    my $seqs = Sequences::FromFasta(shift @$site_files);
    for (@$site_files) {
        my $seq = Sequences::FromFasta($_);
        $seqs = Sequences::Combine($seqs, $seq);
    }   
    Sequences::ToFasta($seqs, $sequences_filename);
        
}
