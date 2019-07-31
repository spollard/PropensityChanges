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
    my $sites_dir = $ARGV[0] // "../data/100seq_test_g100/sites/";

    summarize_acceptabilities($sites_dir);
}

sub summarize_acceptabilities {
    my ($sites_dir) = @_;
    my $acceptabilities = {};

    opendir(my $dir, $sites_dir) or die $!;
    my $site_files = [grep { -f $_} map {"$sites_dir/$_"} grep {/^\d+\.stats\.txt$/} readdir($dir)];
    $site_files = [sort {($a =~ /(\d+)\.stats/)[0] <=> ($b =~ /(\d+)\.stats/)[0]} @$site_files];
    
    # Need to sort the file names by number, not by string
    for my $site_file (@$site_files){
        open my $fh, "<", $site_file or die $!;
        
        # burn header line
        scalar(<$fh>);

        while (my $line = <$fh>){
            chomp $line;
            next if not $line;
            my $fields = [split " ", $line];
            my $node_id = shift @$fields;
            
            my $node_acc = "";
            while (my ($i, $value) = each(@$fields)){
                $node_acc .= ($value >= 0.3) ? $amino_acids->[$i] : '';
            }
            $acceptabilities->{length $node_acc}++;
        }
    }
    
    for my $i (sort keys %$acceptabilities) {
        print "$i: $acceptabilities->{$i}\n";
    }
}
