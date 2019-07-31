# Takes all substitutions from a suboutfile and randomly reassigns them to
# sites and branches that had a substitution. It first reads in all 
# substitutions, then randomly reassigns an end state from the list 
# of previously observed end states for a given start state.  


use strict; 
use warnings;

package Suboutfile;

use Substitutions;

unless(caller) {
	my $suboutfile = $ARGV[0] // "C:/Users/Stephen/Documents/Lab work/Projects/Protein Evolution/Convergence/Data/suboutfile_highprob_ASRV_updatebls7_new";
	my $suboutfile_out = $ARGV[1] // "C:/Users/Stephen/Documents/Lab work/Projects/Protein Evolution/Convergence/Data/suboutfile_highprob_ASRV_updatebls7_new_T";
	my $state = $ARGV[2] // "T";
	
	FilterByAncestralState($suboutfile, $suboutfile_out, $state);
}

sub FilterByAncestralState {
	my ($suboutfile,$suboutfile_out, $state) = @_;
	
	my $substitutions = Substitutions::FromSuboutfile($suboutfile);
	
	for my $subs (values %$substitutions){
		@$subs = grep {substr($_,0,1) eq $state} @$subs;	
	}
	
	Substitutions::ToSuboutfile($substitutions, $suboutfile_out);
}
