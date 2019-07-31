use strict;
use warnings;

use Data::Dumper;

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;


unless (caller){
    my $sites_dir = $ARGV[0] // "../data/simulated_100_sites_100_taxa_g1000_sampled/sites/";
    my $probs_filename = $ARGV[1] // "../data/simulated_100_sites_100_taxa_g1000_sampled/probs";

    collect_acceptability_probabilities($sites_dir, $probs_filename);
}

sub collect_acceptability_probabilities {
    my ($sites_dir, $probs_filename) = @_;
    my $acceptabilities = [];

    opendir(my $dir, $sites_dir) or die $!;
    my $site_files = [grep { -f $_} map {"$sites_dir/$_"} grep {/^\d+\.stats\.txt$/} readdir($dir)];
    $site_files = [sort {($a =~ /(\d+)\.stats/)[0] <=> ($b =~ /(\d+)\.stats/)[0]} @$site_files];

    print(scalar(@$site_files));

    # Need to sort the file names by number, not by string
    for my $site_file (@$site_files){
    print $site_file, "\n";
        open my $fh, "<", $site_file or die $!;
        
        # burn header line
        scalar(<$fh>);

        while (my $line = <$fh>){
            chomp $line;
            next if not $line;
            my $fields = [split " ", $line];
            my $node_id = shift @$fields;
            
            while (my ($i, $value) = each(@$fields)){
                push @$acceptabilities, $value;
            }
        }
    }
    
    open my $fh, ">", $probs_filename or die;
    print $fh join("\n", @$acceptabilities)
}
