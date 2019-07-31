# Turn a suboutfile and treeoutfile into nodepairinfo. 
# The new treeoutfile format does not have the ##..## branch labels in the
# tree. It only has node labels.

use strict;
use warnings;

use Data::Dumper;

use Tree;
use Substitutions;

my $treeoutfile =  $ARGV[0]
  // "../Data/PLEX_data_used_for_paper/treeoutfile_averaged.newick";
my $suboutfile = $ARGV[1]
  // "../Data/PLEX_data_used_for_paper/suboutfile_highprob_2014_10_23_new";
my $nodePairInfo =  $ARGV[2]
  // "../Data/PLEX_data_used_for_paper/NodePairInfo_with_sum_bls_old";

open my $nodePairInfo_out, '>', $nodePairInfo;

my $tree = Tree::FromFile($treeoutfile);

my $substitutions = Substitutions::FromSuboutfile($suboutfile);
Substitutions::SplitSubstitutions($substitutions);
Tree::AttachSubstitutions($tree, $substitutions);

print $nodePairInfo_out
"Left_Node\tRight_Node\tDistance\tSameAnc\tDiffAnc\tConvergence\tDivergence\tCAConvergence\tCADivergence\tDAConvergence\tDADivergence\n";
PrintNodePairInfoAndGetDescendants($tree);

sub PrintNodePairInfoAndGetDescendants {
	my ($tree)  = (@_);

	my $descendants = [];

	if ($tree->{Left} and $tree->{Right}) {

		my $left_descendants =PrintNodePairInfoAndGetDescendants($tree->{Left});
		my  $right_descendants =
		  PrintNodePairInfoAndGetDescendants($tree->{Right});

		push @$left_descendants, $tree->{Left};
		push @$right_descendants, $tree->{Right};

		for my $desc_left (@$left_descendants) {
			for my $desc_right (@$right_descendants) {
				# Skip siblings
				next
				  if $desc_left == $tree->{Left}
					  and $desc_right == $tree->{Right};

				print $nodePairInfo_out $desc_left->{Name};
				print $nodePairInfo_out "\t" . $desc_right->{Name};
				print $nodePairInfo_out "\t"
				  . ($desc_left->{Distance} + $desc_right->{Distance});

				PrintPairInfo($desc_left, $desc_right);
				print $nodePairInfo_out "\n";
			}
		}

		push @$descendants, @$left_descendants;
		push @$descendants, @$right_descendants;

		for my $desc (@$descendants) {
			$desc->{Distance} += $tree->{Distance};
		}
	}

	return $descendants;
}

sub PrintPairInfo {
	my ($node_1, $node_2) = @_;
	my $subs1 = $node_1->{Substitutions};
	my $subs2 = $node_2->{Substitutions};
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
	print $nodePairInfo_out "\t";
	print $nodePairInfo_out join("\t",
		$SameAnc, $DiffAnc, $SameDes, $DiffDes,
		$SA2SD, $SA2DD, $DA2SD, $DA2DD);
}

sub intersection {
	my ($a, $b) = @_;
	my %a = map { $_ => 1} @$a;
	return grep { $a{$_} } @$b;
}
