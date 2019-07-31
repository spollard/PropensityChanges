# Read a tree file and turns it into a single Newick tree

use strict; 
use warnings;

use Treeoutfile;

package Treeoutfile;

sub ToSingleNewickTree { 
	my ($tree_file, $burnin, $generations) = @_;
	my $single_tree_file;
	if (Treeoutfile::ContainsMultipleTrees($tree_file)) {
		print
		  "Averaging the tree file using $generations generations after $burnin 
		 generations of burnin\n";
		$single_tree_file = $tree_file . "_averaged";
		Treeoutfile::AverageWithinFile($tree_file, $burnin, $generations,
			$single_tree_file);
	}
	else {
		print "Tree file contains a single generation\n";
		$single_tree_file = $tree_file;
	}

	my $newick_tree_file;
	if (Treeoutfile::IsOld($single_tree_file)) {
		print "Converting tree file to Newick format\n";
		$newick_tree_file = $single_tree_file . ".newick";
		Treeoutfile::ToNewickTreeFile($single_tree_file,
			$newick_tree_file);
	}
	else {
		print "Tree file is in Newick format\n";
		$newick_tree_file = $single_tree_file;
	}
	return $newick_tree_file;
}

1;