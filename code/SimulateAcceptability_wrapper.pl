# Simulate a covarion model over a tree
# With a fixed amino acid profile.
#
#
#
# Stephen Pollard 2015-6-4

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;

use Options;

#require "MultithreadCommands.pl";
require "SequentialCommands.pl";

unless (caller){

    my $sites_dir = $ARGV[0] // "../data/simulated_100_sites_100_taxa_5/sites/";
    my $number_of_sites = $ARGV[1] // 100;
    mkdir $sites_dir if not -d $sites_dir;

    my $commands = [];

    for my $site (1 .. $number_of_sites) {
        mkdir $sites_dir . $site if not -d $sites_dir . $site;

        my $options = {
            tree_filename   =>   "../data/100.tree",
            tree_out_filename   =>   $sites_dir . $site . "/tree_out",
            sequences_filename    =>  $sites_dir . $site . "/sequences",
            hidden_states_filename  => $sites_dir . $site . "/hidden_states",
            switch_rate  => 0.1,
            sub_rate  =>   0.1,
            q_a       =>   0.1,
            q_u        =>    0.001,
            max_branch_length => 0.05,
        };
        Options::Write($options, $sites_dir . $site . "/ctrl");

        my $command = "perl SimulateSingleAcceptabilityModel.pl $sites_dir$site/ctrl";
        push @$commands, $command;
    }
    
    push @$commands, "perl ConcatenateSites.pl $sites_dir";
    push @$commands, "perl ConcatenateHiddenStates.pl $sites_dir";
    
    #MultithreadCommands($commands);
    SequentialCommands($commands);
}

