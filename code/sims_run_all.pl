# Read a tree file, draw a phylogram using GD::Simple
# color subtree

use strict;
use warnings;

require "GetUniqueDirectory.pl";

unless (caller) {
	my $n = $ARGV[0] // 2;
	my $run_directory = $ARGV[1] // GetUniqueDirectory("../data/");
    my $tree = $ARGV[2] // "../data/treeoutfile_averaged.newick";
    
    mkdir $run_directory if not -d $run_directory;  
    
    my $sites_directory = $run_directory . "sites/";
    
    my $fasta = $run_directory . "all.fasta";
    my $substitutions = $run_directory . "substitutions";
    my $subs_no_multiple = $substitutions . "_no_multiple";
    my $node_pairs = $run_directory . "npi";
    
    my $results = $run_directory;
    $results =~ s/data/results/;
    mkdir $results if not -d $results;
    my $graph = $results . "all.pdf";
    my $same = $results . "same.pdf";
    
    my $log_file = $run_directory . "log.txt";
    
    
	my $commands = ["perl SimulateFixedCovarion_wrapper.pl $sites_directory $n > $log_file",
	"perl ConcatenateSites.pl $sites_directory $fasta >> $log_file",
	"perl SequencesToSuboutfile.pl $fasta $substitutions >> $log_file",
	"perl Suboutfile_9FilterForMultipleSubs.pl $substitutions $sites_directory $subs_no_multiple >> $log_file",
	"perl Suboutfile_ToNodePairInfoRecursively.pm $tree $subs_no_multiple $node_pairs >> $log_file",
    "Rscript PlotNodePairInfo_wrapper.R $node_pairs $graph all 0 >> $log_file",
    "Rscript PlotNodePairInfo_wrapper.R $node_pairs $same same 0 >> $log_file"];
    
    SequentialCommands($commands);
}


sub SequentialCommands {
	my ($commands) = @_;
	
	my $completed_commands = [];
	while (@$commands) {
		my $command = shift @$commands; 
		if (system($command)) {
			print STDERR "This command failed:\n\t$command\n";
			print STDERR "These commands completed successfully:\n\t", join("\n\t", @$completed_commands), "\n";
			print STDERR "These commands were not run:\n", join("\n\t", @$commands), "\n";
			die "Error code: $?";
		}
		push @$completed_commands, $command;
	}

}
