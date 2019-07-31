# Counts the number of substitutions per site

use strict;
use warnings;

use Substitutions;

#my $suboutfile = "../Data/suboutfile_test";
my $suboutfile = "C:/Users/Stephen/Documents/Lab work/Projects/AminoAcidClasses/data/suboutfile_g4000_p90_new";
my $sitesFile = "C:/Users/Stephen/Documents/Lab work/Projects/AminoAcidClasses/data/number_of_subs_per_site";
my $sites = Substitutions::CountNumberOfSubsPerSite($suboutfile);
Substitutions::PrintNumberOfSubsPerSite($sites, $sitesFile);
