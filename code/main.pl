
use Cwd qw( abs_path );
use File::Basename qw( dirname );
use lib dirname(abs_path($0));

use Data::Dumper;

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;

use Options;
use Tree;

my $options = Options::Read("options.txt");


my $tree = Tree::FromFile("../data/100.tree");

#print Tree::ToString($tree);exit;

Tree::recurse_pre($tree, \&dividedistance);

sub dividedistance{
    my ($node) = @_;
    $node->{Distance} = ($node->{Distance}//0) / 10.0;
    print("$node->{Name}, $node->{Distance}","\n");
}

Tree::ToFile($tree, "../data/100.tree.small");
