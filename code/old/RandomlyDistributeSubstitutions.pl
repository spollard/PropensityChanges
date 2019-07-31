# Randomly distributes a set number of substitutions over the tree
#
#
# Stephen Pollard
# 7/24/2015

use strict;
use warnings FATAL => 'all';

#use List::Util 'sum';
#use List::MoreUtils 'any';

#use Data::Dumper;

use Tree;

my $number_of_substitutions = 1;
my $tree_filename = "../data/tree";

srand(0);

my $tree = Tree::FromNewick($tree_filename);

my $length = Tree::CalculateLength($tree);

my $subs = GenerateRandomDistances($length, $number_of_substitutions);

# Must use a hash because I want the set of branches.
my $branches_with_subs = {};
MapSubDistanceToTree($tree, 0, $subs, $branches_with_subs);

my $branch_pair_subs = MakeAllPairs($branches_with_subs);

# Perhaps this should be refactored out, so that it can be used for 
# any tree 
my $distances = FindAllDistances($branch_pair_subs);




sub MakeAllPairs {
	my ($branches_with_subs) = @_;
	
	my $pairs = [];
	for my $i (0 .. $#$branches_with_subs-1) {
		for my $j ($i + 1 .. $#$branches_with_subs){
			push @$pairs, [$branches_with_subs->[$i],$branches_with_subs->[$j]]; 
		}	
	}	
	return $pairs;
}


# flattens the tree in top, left, right order then calculates the upper bound
# of the total length that is assigned to each branch.
sub MapSubDistanceToTree {
	my ($tree, $current_length, $sub_distances, $branches_with_subs) = @_;

	for my $sub_dist (@$sub_distances) {
		if (    $sub_dist >=  $current_length
			and $sub_dist <= $current_length + $tree->{Distance})
		{
			push @{$tree->{sub_dists}}, $sub_dist - $current_length;
			$branches_with_subs->{$tree->{Name}} = 1;
		}
	}

	$current_length += $tree->{Distance};

	if ($tree->{Left} and $tree->{Right}) {
		$current_length = MapSubDistanceToTree(
			$tree->{Left}, $current_length,
			$sub_distances, $branches_with_subs
		);
		$current_length = MapSubDistanceToTree(
			$tree->{Right}, $current_length,
			$sub_distances, $branches_with_subs
		);
	}

	return $current_length;
}

sub GenerateRandomDistances {
	my ($length, $number_of_subs) = @_;

	my $distances = [];

	for (1 .. $number_of_subs) {
		push @$distances, rand($length);
	}

	return $distances;
}

