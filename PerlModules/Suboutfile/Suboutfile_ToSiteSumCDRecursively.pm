# Turn a suboutfile and treeoutfile into site CD.
# The new treeoutfile format does not have the ##..## branch labels in the
# tree. It only has node labels.

use strict;
use warnings;

use Data::Dumper;

use lib "../../../../Perl_modules"; # if run from this directory...
use lib "../../../Perl_modules"; # if used as a module
use Tree;
use Substitutions;

package Suboutfile;

unless (caller) {
	my $SiteCDs = [
		undef,
		{
			"SameAnc" =>1,
			"DiffAnc"=>4,
			"SameDes"=>6,
			"DiffDes"=>7,
			"SA2SD"=>8,
			"SA2DD"=>9,
			"DA2SD"=>10,
			"DA2DD"=>11,
		}
	];
	PrintSiteCD($SiteCDs, 'temp_sitecds');
#	exit;

	my $treeoutfile = $ARGV[0] // "../../Data/treeoutfile_new";

	my $suboutfile = $ARGV[1]
	  // "../../Data/suboutfile_highprob_ASRV_updatebls7_new";
	my $SiteCDs_file = $ARGV[2]// "../../Data/SiteCDs";

	ToSiteSumCD($treeoutfile, $suboutfile, $SiteCDs_file);
}

sub ToSiteSumCD {
	my ($treeoutfile, $suboutfile, $SiteCDs_file) = @_;
	my $tree = Tree::FromFile($treeoutfile);

	my $substitutions = Substitutions::FromSuboutfile($suboutfile);

	Substitutions::SplitSubstitutions($substitutions);
	Tree::AttachSubstitutions($tree, $substitutions);

	my $SiteCDs = [];
	UpdateSiteCDAndGetDescendants($tree, $SiteCDs);

	PrintSiteCD($SiteCDs, $SiteCDs_file);
}

sub UpdateSiteCDAndGetDescendants {
	my ($tree, $SiteCDs)  = (@_);

	my $descendants = [];

	if ($tree->{Left} and $tree->{Right}) {
		$tree->{Left}->{Distance_from_common_ancestor} = 0;
		$tree->{Right}->{Distance_from_common_ancestor} = 0;

		my $left_descendants =
		  UpdateSiteCDAndGetDescendants($tree->{Left}, $SiteCDs);
		my  $right_descendants =
		  UpdateSiteCDAndGetDescendants($tree->{Right}, $SiteCDs);

		push @$left_descendants, $tree->{Left};
		push @$right_descendants, $tree->{Right};

		for my $desc_left (@$left_descendants) {
			for my $desc_right (@$right_descendants) {

				# Skip siblings
				next
				  if $desc_left == $tree->{Left}
					  and $desc_right == $tree->{Right};

				UpdateSiteCD($desc_left, $desc_right, $SiteCDs);
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

sub UpdateSiteCD {
	my ($node_1, $node_2, $SiteCDs) = @_;

	my $subs1 = $node_1->{Substitutions};
	my $subs2 = $node_2->{Substitutions};
	my $commonSites = [intersection([keys %$subs1],	[keys %$subs2])];

	# Sites should be 1 indexed
	foreach my $site (@$commonSites) {
		if ( $subs1->{$site}->{from} eq $subs2->{$site}->{from} ) {
			$SiteCDs->[$site]->{SameAnc}++;
			if ( $subs1->{$site}->{to} eq $subs2->{$site}->{to} ) {
				$SiteCDs->[$site]->{SameDes}++;
				$SiteCDs->[$site]->{SA2SD}++;
			}
			else {
				$SiteCDs->[$site]->{DiffDes}++;
				$SiteCDs->[$site]->{SA2DD}++;
			}
		}
		else {
			$SiteCDs->[$site]->{DiffAnc}++;
			if ( $subs1->{$site}->{to} eq $subs2->{$site}->{to} ) {
				$SiteCDs->[$site]->{SameDes}++;
				$SiteCDs->[$site]->{DA2SD}++;
			}
			else {
				$SiteCDs->[$site]->{DiffDes}++;
				$SiteCDs->[$site]->{DA2DD}++;
			}
		}
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
			print $SiteCDs_out join("\t", $site, map {$_ // 0} @{$SiteCD}{@$keys})."\n";
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
