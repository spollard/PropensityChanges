use strict;
use warnings;

use Data::Dumper;

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;

use Sequences;
use Samplers;

require "amino_acids.pl";
our $amino_acids;

unless (caller){
    my $sites_dir = $ARGV[0] // "../data/simulated_100x100_g10k_sampled/sites/";
    my $fasta_out = $ARGV[1] // "../data/simulated_100x100_g10k_sampled/acceptabilities_3kburnin.fasta";

    collect_probably_acceptabilities($sites_dir, $fasta_out);
}

sub collect_probably_acceptabilities {
    my ($sites_dir, $fasta_out) = @_;
    my $acceptabilities = {};

    opendir(my $dir, $sites_dir) or die $!;
    my $site_files = [grep { -f $_} map {"$sites_dir/$_"} grep {/^\d+\.stats\.txt$/} readdir($dir)];
    $site_files = [sort {($a =~ /(\d+)\.stats/)[0] <=> ($b =~ /(\d+)\.stats/)[0]} @$site_files];

    print(scalar(@$site_files));

    # Need to sort the file names by number, not by string
    for my $site_file (@$site_files){
    print $site_file, "\n";
        open my $fh, "<", $site_file or die $!;
        
        # burn header line
        scalar(<$fh>);

        while (my $line = <$fh>){
            chomp $line;
            next if not $line;
            my $fields = [split " ", $line];
            my $node_id = shift @$fields;
            
            if (exists $acceptabilities->{$node_id}) {
                $acceptabilities->{$node_id}  .= "|";
            }
            
            while (my ($i, $value) = each(@$fields)){
                $acceptabilities->{$node_id} .= ($value >= 0.3) ? $amino_acids->[$i] : '';
            }
        }
    }


    Sequences::ToFasta($acceptabilities, $fasta_out);

}
