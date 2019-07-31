# Tree Functions
# I only use Newick trees and so I don't need Newick in the name

use Tree;

my $filename_in = $ARGV[0] // "AC_cCDS_10heuristic.trimal.3rd.treefix.tre";
my $filename_out = $ARGV[1] // $filename_in .".internal_names";

my $tree = Tree::FromFile($filename_in);
Tree::GenerateNamesForInternalNodes_orderIndependent($tree);

open my $fh_out, ">", $filename_out;

print $fh_out Tree::ToString($tree);

print "\n";
print "@$Tree::names";

