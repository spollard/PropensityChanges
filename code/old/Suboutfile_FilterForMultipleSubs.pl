# THIS FILE WILL NOT WORK IF THE MULTIPLE SUBSTITUTION FILE NAMING SCHEME
# IS CHANGED!

# CollectMultipleSubstitutions() must change if the multiple substitution file 
# name or format changes

 
use strict;
use warnings;

use Data::Dumper;

use Substitutions;

unless (caller) {
	my $substitutions_filename = $ARGV[0] // "../data/2015-08-13T12_03_44_0/substitutions";
	my $sites_dir = $ARGV[1] // "../data/2015-08-13T12_03_44_0/" ;
	my $substitutions_out_filename = $ARGV[2] //"../data/2015-08-13T12_03_44_0/substitutions_no_multiple_subs";
	
	FilterForMultipleSubs($substitutions_filename, $sites_dir, $substitutions_out_filename);
	
}

sub FilterForMultipleSubs {
	my ($substitutions_filename, $sites_dir, $substitutions_out_filename) = @_;

	my $multiple_subs = CollectMultipleSubstitutions($sites_dir);
	my $substitutions = Substitutions::FromSuboutfile($substitutions_filename);
	Substitutions::RemoveMultipleSubs($substitutions, $multiple_subs);
	Substitutions::ToSuboutfile($substitutions, $substitutions_out_filename);
}


sub CollectMultipleSubstitutions {
	my ($sites_dir) = @_;
	
	my $multiple_subs = {}; #{branch_name =>[1, 2, 5]}
	
	my $multiple_subs_filenames = [<$sites_dir*.multiple_subs>];
	
	for my $file (@$multiple_subs_filenames) {
		$file =~ /site(\d+)_/;
		my $site = $1;
		open my $f, "<", $file;
		while(<$f>){
			chomp;
			my $name = $_;
			if (not defined $multiple_subs->{$name}){
					$multiple_subs->{$name} = [];
			}
			push @{$multiple_subs->{$name}}, $site;
		}
	}
	
	return $multiple_subs;
}

package Substitutions;
sub RemoveMultipleSubs {
	my ($substitutions, $multiple_subs) = @_;
		
	# substitutions are currently flat
	
	Split($substitutions);
	
	while (my ($name, $subs) = each %$substitutions) {
		next if not defined $multiple_subs->{$name};
		for my $site (@{$multiple_subs->{$name}}){
			print "Removing sub from branch $name at site $site\n";
			delete $subs->{$site} if defined $subs->{$site};
		}
	}
	
	# Make them flat again
	Join($substitutions);
}
