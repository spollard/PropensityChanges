use strict;
use warnings;

use Data::Dumper;

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;

use Sequences;
use Samplers;

require "amino_acids.pl";
our $amino_acids;

unless (caller){
    my $sites_dir = $ARGV[0] // "../data/100seq_test_g100/sites/";
    my $fasta_out = $ARGV[1] // "../data/100seq_test_g100/simulated.fasta";

    simulate_from_acceptabilities($sites_dir, $fasta_out);
}

sub simulate_from_acceptabilities {
    my ($sites_dir, $fasta_out) = @_;
    my $sequences = {};

    opendir(my $dir, $sites_dir) or die $!;
    my $site_files = [grep { -f $_} map {"$sites_dir/$_"} grep {/^\d+\.stats\.txt$/} readdir($dir)];
    $site_files = [sort {($a =~ /(\d+)\.stats/)[0] <=> ($b =~ /(\d+)\.stats/)[0]} @$site_files];
    
    # Need to sort the file names by number, not by string
    for my $site_file (@$site_files){
        open my $fh, "<", $site_file or die $!;
        
        # burn header line
        scalar(<$fh>);

        while (my $line = <$fh>){
            chomp $line;
            next if not $line;
            my $fields = [split " ", $line];
            my $node_id = shift @$fields;
            
            $sequences->{$node_id} .= $amino_acids->[Samplers::weighted_random_choice($fields)]
        }
    }


    Sequences::ToFasta($sequences, $fasta_out);

}
