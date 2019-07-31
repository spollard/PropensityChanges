# Turn a directory of suboutfiles and treeoutfile into a directory of
# nodepairinfos.

use strict;
use warnings;

use Suboutfile_OldToNew;

unless (caller) {
	my $treeoutfile = $ARGV[0] // "treeoutfile_averaged";
	my $suboutfile_dir = $ARGV[1]// "suboutfile_splits/";
	my $nodePairInfo_dir = $ARGV[2]// "NodePairInfos/";
	my $nodes_file = $ARGV[3] // "nodes";

	mkdir $nodePairInfo_dir if not -d $nodePairInfo_dir;

	my $suboutfiles = [<$suboutfile_dir*>];

	# All the suboutfiles need to be updated before they can be used
	for my $suboutfile (@$suboutfiles) {
		my $suboutfile_new = $suboutfile . "_new";
		Suboutfile::OldToNew($suboutfile, $nodes_file, $suboutfile_new);
	}
	exit;

	# The above loop must end before the next loop can begin.
	for my $s ((1120*2) .. $#$suboutfiles) {
		my $suboutfile = $suboutfiles->[$s];
		next if $suboutfile =~ /_new/;
		my $nodepairinfo_file = $suboutfile;
		$nodepairinfo_file =~ s/$suboutfile_dir/$nodePairInfo_dir/;
		$suboutfile .= "_new";
		print "working on file $suboutfile\n";
		if ($s % 18 == 0) {
			system
"perl Suboutfile_ToNodePairInfoRecursively.pm $treeoutfile $suboutfile $nodepairinfo_file";
		}
		else {
			system
"perl Suboutfile_ToNodePairInfoRecursively.pm $treeoutfile $suboutfile $nodepairinfo_file &";
		}
	}
}
