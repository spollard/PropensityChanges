# Turn a suboutfile and treeoutfile into nodepairinfo.
# The new treeoutfile format does not have the ##..## branch labels in the
# tree. It only has node labels.

# Cragganmore cannot hold all of this information in memory at once. 
# I think it is the substitiutions 


use strict;
use warnings;

use Data::Dumper;

# There is another 'Tree.pm'
use lib ".";
use Tree;
use Substitutions;

package Suboutfile;

unless (caller) {

	#my $treeoutfile = $ARGV[0]
	#  // "../Data/PLEX_data_used_for_paper/treeoutfile_averaged.newick";

	my $treeoutfile = $ARGV[0] // "treeoutfile_new";

	my $suboutfile = $ARGV[1]// "suboutfile";
	my $burnin = $ARGV[2] // 1000;
	my $gens = $ARGV[3] // 4000;

	my $nodePairInfo_file = $ARGV[4]// "NodePairInfo_manyGens";

	ToNodePairInfo($treeoutfile, $suboutfile, $burnin, $gens,
		$nodePairInfo_file);
}

sub ToNodePairInfo {
	my ($treeoutfile, $suboutfile, $burnin, $gens, $nodePairInfo_file) = @_;

	my $tree = Tree::FromFile($treeoutfile);

	my $substitutions =
	  Substitutions::ManyGensFromSuboutfile($suboutfile, $burnin, $gens);
	# This blows up in memory
#	for my $subs (values $substitutions) {
#		Substitutions::SplitSubstitutions($subs);	
#	}
	# This will work for many generations of substitutions also
	Tree::AttachSubstitutions($tree, $substitutions);

	open my $nodePairInfo_out, '>', $nodePairInfo_file;

	print $nodePairInfo_out join("\t",
		"Left_Node", "Right_Node", "Distance","Sum_Branch_Lengths",
		map {"Gen_" . $_} $burnin..($burnin+$gens-1))
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

# This needs to be changed to handle multiple generations of substitutions
sub PrintPairInfo {
	my ($node_1, $node_2, $nodePairInfo_out) = @_;

	my $distance_from_eachother =
	  $node_1->{Distance_from_common_ancestor} +
	  $node_2->{Distance_from_common_ancestor};
	my $sum_branch_lengths = $node_1->{Distance} + $node_2->{Distance};
	print $nodePairInfo_out join("\t",
		$node_1->{Name}, $node_2->{Name},
		$distance_from_eachother, $sum_branch_lengths,
	);
	for my $gen (0 .. $#{$node_1->{Substitutions}}){
		my $subs1 = $node_1->{Substitutions}->[$gen];
		my $subs2 = $node_2->{Substitutions}->[$gen];
		my $commonSites = [intersection([keys %$subs1],	[keys %$subs2])];

		my $SameAnc = 0;
		my $DiffAnc = 0;

		my $SameDes = 0; #conv
		my $DiffDes = 0; #div

		my $SA2SD =
		  0; # Same Ancestor to Same Descendant, Common ancestor convergence
		my $SA2DD = 0; # Common ancestor divergence
		my $DA2SD = 0; # Different Ancestor to Same Descendent
		my $DA2DD = 0; # Different Ancestor to Different Descendent

		foreach my $site (@$commonSites) {
			if ( $subs1->{$site}->{from} eq $subs2->{$site}->{from} ) {
				$SameAnc++;
				if ( $subs1->{$site}->{to} eq $subs2->{$site}->{to} ) {
					$SameDes++;
					$SA2SD++;
				}
				else {
					$DiffDes++;
					$SA2DD++;
				}
			}
			else {
				$DiffAnc++;
				if ( $subs1->{$site}->{to} eq $subs2->{$site}->{to} ) {
					$SameDes++;
					$DA2SD++;
				}
				else {
					$DiffDes++;
					$DA2DD++;
				}
			}
		}

		print $nodePairInfo_out "\t"
		  . join(",",
			$SameAnc,$DiffAnc,$SameDes,$DiffDes,$SA2SD,$SA2DD,$DA2SD,$DA2DD)
		  ;
	}
	print $nodePairInfo_out "\n";
}

sub intersection {
	my ($a, $b) = @_;
	my %a = map { $_ => 1} @$a;
	return grep { $a{$_} } @$b;
}

