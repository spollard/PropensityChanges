# Simply pass all the commands you want to call in order.
#
# The SequentialCommands subroutine will block until all commands have 
# finished running. 
# 
# Only tested on Windows 10 and Ubuntu 12.10 and 16.04 
#
# Stephen Pollard 12/1/2016
# 
# Usage: 
# require "SequentialCommands.pl"
# 
# # $commands is a ref to an array of command strings 
#
# SequentialCommands($commands); 

use strict;
use warnings;
use diagnostics;

unless(caller) {
	# First generate all the commands to be run
	# This for loop could loop over a parameter changing
	my $commands = [];
	
	for (0 .. 10) {
		push @$commands, "echo $_";
		if ($_ == 7) {push @$commands, "intentional_fail_command";}
	}
	
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
			print STDERR "These commands were not run:\n\t", join("\n\t", @$commands), "\n";
			die "Error code: $?";
		}
		push @$completed_commands, $command;
	}
}