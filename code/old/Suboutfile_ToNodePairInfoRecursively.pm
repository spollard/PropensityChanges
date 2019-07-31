# Turn a suboutfile and treeoutfile into nodepairinfo.
# The new treeoutfile format does not have the ##..## branch labels in the
# tree. It only has node labels.

use strict;
use warnings;

use Data::Dumper;

use Tree;
use Substitutions;

package Suboutfile;

unless (caller) {

	#my $treeoutfile = $ARGV[0]
	#  // "../Data/PLEX_data_used_for_paper/treeoutfile_averaged.newick";

	my $treeoutfile = $ARGV[0] // "../../Data/treeoutfile_new";

	my $suboutfile = $ARGV[1]
	  // "../../Data/suboutfile_highprob_ASRV_updatebls7_new";
	my $nodePairInfo_file = $ARGV[2]// "../../Data/NodePairInfo_CO1_only";

	ToNodePairInfo($treeoutfile, $suboutfile, $nodePairInfo_file);
}

sub ToNodePairInfo {
	my ($treeoutfile, $suboutfile, $nodePairInfo_file) = @_;
	my $tree = Tree::FromFile($treeoutfile);
    
	my $substitutions = Substitutions::FromSuboutfile($suboutfile);

	Substitutions::SplitSubstitutions($substitutions);
	Tree::AttachSubstitutions($tree, $substitutions);

	open my $nodePairInfo_out, '>', $nodePairInfo_file;

	print $nodePairInfo_out join("\t",
		"Left_Node", "Right_Node", "Distance",
		"Sum_Branch_Lengths", "SameAnc","DiffAnc",
		"Convergence", "Divergence", "CAConvergence",
		"CADivergence","DAConvergence", "DADivergence")
	  . "\n";
	PrintNodePairInfoAndGetDescendants($tree, $nodePairInfo_out);

}

sub PrintNodePairInfoAndGetDescendants {
	my ($tree, $nodePairInfo_out)  = (@_);

	my $descendants = [];

	if ($tree->{Left} and $tree->{Right}) {
		$tree->{Left}->{Distance_from_common_ancestor} = 0;
		$tree->{Right}->{Distance_from_common_ancestor} = 0;

		my $left_descendants =
		  PrintNodePairInfoAndGetDescendants($tree->{Left}, $nodePairInfo_out);
		my  $right_descendants =
		  PrintNodePairInfoAndGetDescendants($tree->{Right}, $nodePairInfo_out);

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

	my $subs1 = $node_1->{Substitutions};
	my $subs2 = $node_2->{Substitutions};
	my $commonSites = [intersection([keys %$subs1],	[keys %$subs2])];

	my $sameAncestralResidues = 0;
	my $differentAncestralResidues = 0;

	my $sameDescendantResidues = 0;
	my $differentDescendantResidues = 0;

	my $sameAncestralResiduesAndSameDescendantResidues =0;
	my $sameAncestralResiduesAndDifferentDescendantResidues =0;
	my $differentAncestralResiduesAndSameDescendantResidues =0;
	my $differentAncestralResiduesAndDifferentDescendantResidues =0;

	foreach my $site (@$commonSites) {
		if ( $subs1->{$site}->{from} eq $subs2->{$site}->{from} ) {
			$sameAncestralResidues++;
			if ( $subs1->{$site}->{to} eq $subs2->{$site}->{to} ) {
				$sameDescendantResidues++;
				$sameAncestralResiduesAndSameDescendantResidues++;
			}
			else {
				$differentDescendantResidues++;
				$sameAncestralResiduesAndDifferentDescendantResidues++;
			}
		}
		else {
			$differentAncestralResidues++;
			if ( $subs1->{$site}->{to} eq $subs2->{$site}->{to} ) {
				$sameDescendantResidues++;
				$differentAncestralResiduesAndSameDescendantResidues++;
			}
			else {
				$differentDescendantResidues++;
				$differentAncestralResiduesAndDifferentDescendantResidues++;
			}
		}
	}

	my $distance_from_eachother =
	  $node_1->{Distance_from_common_ancestor} +
	  $node_2->{Distance_from_common_ancestor};
	my $sum_branch_lengths = $node_1->{Distance} + $node_2->{Distance};
	print $nodePairInfo_out join("\t",
		$node_1->{Name},
		$node_2->{Name},
		$distance_from_eachother,
		$sum_branch_lengths,
		$sameAncestralResidues,
		$differentAncestralResidues,
		$sameDescendantResidues,
		$differentDescendantResidues,
		$sameAncestralResiduesAndSameDescendantResidues,
		$sameAncestralResiduesAndDifferentDescendantResidues,
		$differentAncestralResiduesAndSameDescendantResidues,
		$differentAncestralResiduesAndDifferentDescendantResidues)
	  . "\n";
}

sub intersection {
	my ($a, $b) = @_;
	my %a = map { $_ => 1} @$a;
	return grep { $a{$_} } @$b;
}

1;
