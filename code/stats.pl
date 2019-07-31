

use strict;
use warnings FATAL=>'all';
use diagnostics;

use Data::Dumper;
#~ use Set::Tiny;


our $PerlModulesDir;
BEGIN { require "./PerlModulesDir.pl"; }
use lib $PerlModulesDir;

use Options;
use Statistics qw(average);

require "./amino_acids.pl";
our $amino_acids;

unless(caller){
    # Default options
    my $options = Options::Read("stats.options");
    
    # Provide options files on the command line. Later files overwrite newer files
    for my $file (@ARGV){
        print "Reading options from $file\n";
        %$options = (%$options, %{Options::Read($file)});
    }
    
    AcceptableMCMC_stats($options);
}


sub AcceptableMCMC_stats {
    my ($options) = @_;
    
    my $model_filename = $options->{out_directory}."model";
    my $sites_dir = $options->{out_directory}."sites/";
    my $burnin = $options->{burnin} // 100;
    my $skipped_columns = 0;
    my $summary_out_filename = $options->{out_directory} . "summary.txt";
    
    opendir(my $dir, $sites_dir) or die $!;
    my $site_files = [grep { -f $_} map {"$sites_dir/$_"} grep {/^\d+$/} readdir($dir)];
    
    my $acceptable_distribution = [];
    
    for my $site_file (@$site_files){
        open my $fh, "<", $site_file or die $!;
        
        # Get node ids
        my $line = scalar(<$fh>);
        chomp $line;
        my $node_ids = [ split("\t", $line)];
        $node_ids = [@{$node_ids}[$skipped_columns..$#$node_ids]];
        
        
        
        my $node_acceptables = {};
        my $generations = 0;
        
        
        while (my $line = <$fh>){
            next if $. < $burnin;
            chomp $line;
            my $parts = [split("\t", $line)];
            $parts = [@{$parts}[$skipped_columns..$#$parts]];
            
            for my $i (0..$#$parts){
                my $acceptables = [$parts->[$i] =~ /(\w)/g];
                
                $acceptable_distribution->[scalar(@$acceptables)]++;
                
                for my $a (@$acceptables){
                    $node_acceptables->{$node_ids->[$i]}->{$a}++;
                }
            }
            $generations++;
        }
        
        print "Read $generations generations\n";
        
        
        open my $fh_out, ">", $site_file . ".stats.txt";
        print $fh_out join("\t", "Node_id", @$amino_acids) . "\n";
        for my $id (@$node_ids) {
            print($fh_out join("\t", $id,map {($_//0)/$generations} @{$node_acceptables->{$id}}{@$amino_acids}) . "\n")
        }
        
    }

    
        
    open my $fh_out_summary, ">", $summary_out_filename;
    print $fh_out_summary join("\t", keys @$acceptable_distribution) . "\n";
    print $fh_out_summary join("\t", map({$_//0} @$acceptable_distribution)) . "\n";
    
    
    open my $fh, "<", $model_filename or die $!;
    print $fh_out_summary scalar(<$fh>);
    
    my $parameters = [];
    while (my $line = <$fh>){
        chomp $line;
        next if not $line;
        my $fields = [split " ", $line];
        while (my ($i, $value) = each(@$fields)){
            push(@{$parameters->[$i]}, $value);
        }
    }
    print $fh_out_summary join("\t", map {average($_)} @$parameters),"\n";
}

