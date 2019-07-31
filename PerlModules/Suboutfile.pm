# Convert an old suboutfile with branch labels (e.g. ##..##) to new suboutfiles
# with node labels (e.g. Node_136 or Canis_lupus_8234).

use strict;
use warnings;

package Suboutfile;


use File::Basename;
use lib dirname(__FILE__) . "/Suboutfile/";
use Suboutfile;
use Suboutfile_FilterBySingleSubstitutionPosteriorProbability;
use Suboutfile_OldToNew;
use Suboutfile_ToNodePairInfoRecursively;

unless (caller) {
	exit;
	print IsOld("../Data/pipeline/suboutfile_highprob_ASRV_updatebls7_new_duplicated");
	print ContainsMultipleGenerations("../Data/pipeline/suboutfile_highprob_ASRV_updatebls7_new_duplicated_highprob");
	print ContainsMultipleGenerations("temp");
	exit;
	
	
	my $suboutfile_old = $ARGV[0] // "../suboutfile_highprob_ASRV_updatebls7";
	my $nodes_file = $ARGV[1] // "nodes";
	my $suboutfile_new = $ARGV[2] // ($suboutfile_old . "_new");

	OldToNew($suboutfile_old, $nodes_file, $suboutfile_new);
}
