# Turn a suboutfile and treeoutfile into site CD.
# The new treeoutfile format does not have the ##..## branch labels in the
# tree. It only has node labels.

use strict;
use warnings FATAL => 'all';

use Data::Dumper;

use Tree;
use Substitutions;

package Suboutfile;

unless (caller) {
	my $treeoutfile = $ARGV[0] // "../../Data/treeoutfile_new";

	my $suboutfile = $ARGV[1]
	  // "../../Data/suboutfile_highprob_ASRV_updatebls7_new";
	my $SiteCDs_directory = $ARGV[2] // "../../Data/SiteIndividualCDs/";
	mkdir $SiteCDs_directory if not -d $SiteCDs_directory;
	
	ToSiteIndividualCD($treeoutfile, $suboutfile, $SiteCDs_directory);
}

sub ToSiteIndividualCD {
	my ($treeoutfile, $suboutfile, $SiteCDs_directory) = @_;

	my $tree = Tree::FromFile($treeoutfile);

	my $substitutions = Substitutions::FromSuboutfile($suboutfile);

	Substitutions::SplitSubstitutions($substitutions);
	Tree::AttachSubstitutions($tree, $substitutions);

	PrintSiteIndividualCDAndGetDescendants($tree, $SiteCDs_directory);
}
 
sub PrintSiteIndividualCDAndGetDescendants {
	my ($tree, $SiteCDs_directory)  = (@_);
	
	my $descendants = [];

	if ($tree->{Left} and $tree->{Right}) {
		$tree->{Left}->{Distance_from_common_ancestor} = 0;
		$tree->{Right}->{Distance_from_common_ancestor} = 0;

		my $left_descendants =
		  PrintSiteIndividualCDAndGetDescendants($tree->{Left}, $SiteCDs_directory);
		my  $right_descendants =
		  PrintSiteIndividualCDAndGetDescendants($tree->{Right}, $SiteCDs_directory);

		push @$left_descendants, $tree->{Left};
		push @$right_descendants, $tree->{Right};

		for my $desc_left (@$left_descendants) {
			for my $desc_right (@$right_descendants) {

				# Skip siblings
				next
				  if $desc_left == $tree->{Left}
					  and $desc_right == $tree->{Right};

				PrintSiteIndividualCD($desc_left, $desc_right,
					$SiteCDs_directory);
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

sub PrintSiteIndividualCD {
	my ($node_1, $node_2, $SiteCDs_directory) = @_;

	my $subs1 = $node_1->{Substitutions};
	my $subs2 = $node_2->{Substitutions};
	
	my $distance_from_eachother =
	  $node_1->{Distance_from_common_ancestor} +
	  $node_2->{Distance_from_common_ancestor};
	my $sum_branch_lengths = $node_1->{Distance} + $node_2->{Distance};
	
	my $commonSites = [intersection([keys %$subs1],	[keys %$subs2])];
	# Sites should be 1 indexed
	foreach my $site (@$commonSites) {

		my $site_out;
		if (not -f $SiteCDs_directory . $site) {
			open $site_out, ">", $SiteCDs_directory . $site;
			print $site_out join("\t",
				"Left_Node", "Right_Node", "Distance",
				"Sum_Branch_Lengths", "SameAnc","DiffAnc",
				"Convergence", "Divergence", "CAConvergence",
				"CADivergence","DAConvergence", "DADivergence")
			  . "\n";

		}
		else {
			open $site_out, ">>", $SiteCDs_directory . $site;
		}

		my $SameAnc = 0;
		my $DiffAnc = 0;

		my $SameDes = 0; #conv
		my $DiffDes = 0; #div

		my $SA2SD =
		  0; # Same Ancestor to Same Descendant, Common ancestor convergence
		my $SA2DD = 0; # Common ancestor divergence
		my $DA2SD = 0; # Different Ancestor to Same Descendent
		my $DA2DD = 0; # Different Ancestor to Different Descendent
		
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

		print $site_out join("\t",
			$node_1->{Name}, $node_2->{Name},
			$distance_from_eachother, $sum_branch_lengths,
			$SameAnc, $DiffAnc,
			$SameDes, $DiffDes,
			$SA2SD, $SA2DD,
			$DA2SD, $DA2DD)
		  . "\n";
	}
}

sub PrintSiteCD {
	my ($SiteCDs, $SiteCDs_file) = @_;

	open my $SiteCDs_out, '>', $SiteCDs_file;
	print $SiteCDs_out join("\t",
		"Site", "SameAnc", "DiffAnc",
		"Convergence", "Divergence","CAConvergence",
		"CADivergence","DAConvergence", "DADivergence")
	  . "\n";

	my $keys = [
		"SameAnc", "DiffAnc","SameDes", "DiffDes",
		"SA2SD", "SA2DD","DA2SD", "DA2DD",
	];
	while (my ($site, $SiteCD) = each(@$SiteCDs)) {

		# $SiteCD will be undefined for sites without substitutions
		if (defined $SiteCD) {
			print $SiteCDs_out
			  join("\t", $site, map {$_ // 0} @{$SiteCD}{@$keys})."\n";
		}
		else {
			print $SiteCDs_out join("\t", $site, (0)x 8)."\n";
		}
	}
}

sub intersection {
	my ($a, $b) = @_;
	my %a = map { $_ => 1} @$a;
	return grep { $a{$_} } @$b;
}

1;
