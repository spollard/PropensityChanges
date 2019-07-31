# Find the average tree (same topology) in a treeoutfile after burnin

use strict;
use warnings;

package Treeoutfile;


unless (caller) {
	my $multiple_tree_file = $ARGV[0] // "../Data/pipeline/treeoutfile_head_new_labels";
	my $burnin = $ARGV[1] // 0;
	my $generations = $ARGV[2] // "inf";
	my $average_tree_file = $ARGV[3] // $multiple_tree_file . "_averaged";
	
	AverageWithinFile($multiple_tree_file, $burnin, $generations, $average_tree_file);
}

sub AverageWithinFile {
	my ($multiple_tree_file, $burnin, $generations, $average_tree_file) = @_;
	
	my $branch_lengths_file = $multiple_tree_file . "_branch_lengths";
	
	ToBranchLengths($multiple_tree_file, $branch_lengths_file);
	
	my $average_branch_lengths = AverageBranchLengths($branch_lengths_file, $burnin, $generations);
	
	PrintAverageTree($multiple_tree_file, $average_branch_lengths, $average_tree_file);
}

sub IsOld {
	my ($tree_file) = @_;
	open my $tree_in, "<", $tree_file or die $!;
	my $line = scalar <$tree_in>;
	my $is_old = ($line =~ /^tree gen/ or $line =~ /\.\./);
	return $is_old;
}

sub ContainsMultipleTrees {
	my ($tree_file) = @_;
	open my $tree_in, "<", $tree_file or die $!;
	my $contains_multiple_trees = 0;
	while (<$tree_in>) {
		if ($. > 1) {
			$contains_multiple_trees = 1;
			last;	
		} 	
	}
	return $contains_multiple_trees;
}


sub AverageBranchLengths {
	my ($branch_lengths_file, $burnin, $generations) = @_;
	
	open my $branch_lengths_in, "<", $branch_lengths_file or die $!;
	open my $branch_lengths_out, ">", $branch_lengths_file . "_averaged" or die $!;
	
	my $header_line = scalar <$branch_lengths_in>;
	my $names = [split "\t", $header_line];
	print $branch_lengths_out $header_line; #keep the header
	
	my $average_branch_lengths = [];
	my $generations_used = 0;
	$. = 0; # Reset the lines read counter
	while (<$branch_lengths_in>) {
		next if ($. <= $burnin);
		last if ($. > $burnin + $generations);

		my $line_branch_lengths = [split];
		
		for (0 .. $#$line_branch_lengths) {
			$average_branch_lengths->[$_] += $line_branch_lengths->[$_];
		}
		$generations_used++;
	}
	
	for (@$average_branch_lengths) {
		$_ /= $generations_used;
		# Perhaps use sprintf here instead
		$_ = sprintf("%05.3",$_); # To prevent scientific notation while
									# printing
	}
	
	
	print $branch_lengths_out join("\t", @$average_branch_lengths);
	
	my %average_branch_lengths;
	@average_branch_lengths{@$names} = @$average_branch_lengths;
	
	# Return the hash ref and not the array ref
	$average_branch_lengths = \%average_branch_lengths;
	return $average_branch_lengths;
}


# The next step is to replace the old distances in the tree with the average 
# distances. This really is a string replacement problem. 
sub PrintAverageTree {
	my ($multiple_tree_file, $average_branch_lengths, $average_tree_file) = @_;
	open my $tree_in, "<", $multiple_tree_file or die $!;
	my $line = scalar <$tree_in>;
	$line =~ s/^tree gen\d* = //;
	# Grab the name and replace the incorrect branch length with the average
	# branch length
	$line =~ s/([^(),]*):\d*\.\d*/$1:$average_branch_lengths->{$1}/g;
	
	open my $tree_out, ">", $average_tree_file;
	
	print $tree_out $line;
}

sub ToNewickTreeFile {
	my ($plex_tree_file, $newick_tree_file) = @_;
	open my $tree_in, "<", $plex_tree_file or die $!;
	open my $tree_out, ">", $newick_tree_file;
	while (my $line = <$tree_in>) {
		$line =~ s/^tree gen\d* = //;
		
		my $newick_tree = PlexTreeToNewickTree($line);
		print $tree_out $newick_tree;
	}
}

sub PlexTreeToNewickTree {
	my ($plex_tree) = @_;
		
	$plex_tree =~ s/\) \d+\.\./\)Node_/g; # This only affects internal branches
	$plex_tree =~ s/ \d+\.\.\d+//g; # Remove label from all branches
	
	$plex_tree =~ s/\);/\)Node_0:0;/; # Append the root label and distance
	my $newick_tree = $plex_tree;
	
	return $newick_tree;
}

sub ToBranchLengths {
	my ($treeoutfile, $distances_file) = @_;
	
	open my $tree_in, $treeoutfile or die $!;
	open my $distances_out, ">", $distances_file;
	my $first = 1;
	while (<$tree_in>) {
		s/^tree gen\d* = //;
		my @matches = ($_ =~ /([^(),]*:[^,();]*)/g);
	
		if ($first) {
			$first = 0;
			for (@matches) {
				/(.*):/;
				print $distances_out "$1\t";
			}
			print $distances_out "\n";
		}
		for (@matches) {
			/:(.*)/;
			print $distances_out "$1\t";
		}
		print $distances_out "\n";
	}
}


1;