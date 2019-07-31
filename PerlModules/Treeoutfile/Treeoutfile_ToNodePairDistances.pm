# Calculates the distance between every pair of branches

use strict;
use warnings;

use Data::Dumper;

use Tree;
use Substitutions;

package Treeoutfile;

unless (caller) {

	my $treeoutfile = $ARGV[0] // "C:/Users/Stephen/Documents/Lab work/Projects/AminoAcidClasses/data/treeoutfile_averaged.newick";

	my $nodePairDistances_filename = $ARGV[1]
	  // "C:/Users/Stephen/Documents/Lab work/Projects/AminoAcidClasses/data/NodePairDistances";

	ToNodePairDistances($treeoutfile, $nodePairDistances_filename);
}

sub ToNodePairDistances {
	my ($treeoutfile, $nodePairDistances_filename) = @_;
	my $tree = Tree::FromFile($treeoutfile);

	
	open my $nodePairDistances_out, '>', $nodePairDistances_filename;

	print $nodePairDistances_out join("\t",
		"Left_Node", "Right_Node", "Distance")
	  . "\n";
	PrintNodePairDistanceAndGetDescendants($tree, $nodePairDistances_out);
	
}

sub PrintNodePairDistanceAndGetDescendants {
	my ($tree, $nodePairInfo_out)  = (@_);

	my $descendants = [];

	if ($tree->{Left} and $tree->{Right}) {
		$tree->{Left}->{Distance_from_common_ancestor} = 0;
		$tree->{Right}->{Distance_from_common_ancestor} = 0;

		my $left_descendants = PrintNodePairDistanceAndGetDescendants($tree->{Left}, $nodePairInfo_out);
		my  $right_descendants =
		  PrintNodePairDistanceAndGetDescendants($tree->{Right}, $nodePairInfo_out);

		push @$left_descendants, $tree->{Left};
		push @$right_descendants, $tree->{Right};

		for my $desc_left (@$left_descendants) {
			for my $desc_right (@$right_descendants) {

				# Skip siblings
				next
				  if $desc_left == $tree->{Left}
					  and $desc_right == $tree->{Right};

				PrintPairInfo($desc_left, $desc_right, $nodePairInfo_out);
			}
		}

		push @$descendants, @$left_descendants;
		push @$descendants, @$right_descendants;

		for my $desc (@$descendants) {
			$desc->{Distance_from_common_ancestor} += $tree->{Distance};
		}
	}

	return $descendants;
}

sub PrintPairInfo {
	my ($node_1, $node_2, $nodePairInfo_out) = @_;

	my $distance_from_eachother =
	  $node_1->{Distance_from_common_ancestor} +
	  $node_2->{Distance_from_common_ancestor};
	print $nodePairInfo_out join("\t",
		$node_1->{Name}, $node_2->{Name},
		$distance_from_eachother)
	  . "\n";
}

1;
