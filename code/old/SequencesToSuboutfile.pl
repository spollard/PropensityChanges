# Read a tree file, draw a phylogram using GD::Simple
# color subtree

use strict;
use warnings;


use Sequences;
use Tree;
use Substitutions;

unless (caller) {
	my $data_dir = "../data/";

	my $sequences_filename = $ARGV[0] // $data_dir . "sites/all.fasta";
	my $tree_filename = $data_dir . "treeoutfile_averaged.newick";
	my $suboutfile_name = $ARGV[1] // $data_dir . "suboutfile";
	
	my $seqs = Sequences::FromFasta($sequences_filename);
	
	my $tree = Tree::FromFile($tree_filename);
	
	Tree::AttachSequences($tree, $seqs);
	Tree::SequencesToSubstitutions($tree);
	
	my $subs = Tree::ToSubstitutions($tree);
		
	Substitutions::ToSuboutfile($subs, $suboutfile_name);
		
}
