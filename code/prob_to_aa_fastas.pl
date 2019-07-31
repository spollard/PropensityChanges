use strict;
use warnings;

use Data::Dumper;

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;

use Sequences;

require "amino_acids.pl";
our $amino_acids;

unless (caller){
    my $sites_dir = $ARGV[0] // "../data/acceptable/sites/";
    my $aa_fastas_dir = $ARGV[1] // "../data/acceptable/";

    prob_to_aa_fastas($sites_dir, $aa_fastas_dir);
}

sub prob_to_aa_fastas {
    my ($sites_dir, $aa_fastas_dir) = @_;
    my $aa_seqs = {};

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
            while (my ($i, $value) = each(@$fields)){
                $aa_seqs->{$amino_acids->[$i]}->{$node_id} .= ($value >= 0.1) ? $amino_acids->[$i] : '0';
            }
        }
    }

    for my $aa (@$amino_acids){
        Sequences::ToFasta($aa_seqs->{$aa}, $aa_fastas_dir . $aa . ".fasta");
    }
}