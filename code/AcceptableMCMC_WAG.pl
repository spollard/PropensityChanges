# Takes a tree and a site and estimates the likelihood of an 
# acceptability model.
# 

# Stephen Pollard
# 2015-6-13

use strict;
use warnings FATAL => 'all';

use File::Copy;

use List::Util qw(sum shuffle);
use List::MoreUtils 'any';

use Set::Tiny;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

use Clone qw(clone);

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;

use Tree;
use Sequences;
use Propensities;
use Samplers;
use Options;


require "MCMC.pl";

require "amino_acids_WAG.pl";
our $amino_acids;

# The problem with 'require'ing the options is that require happens
# at compile time. So you must "do" or "eval" with all the security
#  problems of having config as a scripting language
#~ require "acceptable.options.pl";
#~ our $options;

# Globals for now
my ($ancestral_node_names, $all_node_names);
my $amino_acids_considered = $amino_acids;
my $sampler_stats = [];
my $options;
my $number_of_sites;
my $number_of_branches;
my $number_of_amino_acids = scalar(@$amino_acids_considered);
my $empty_aggregate;

unless(caller){
    
    # Default options
    $options = Options::Read("AcceptableMCMC_WAG.options");

    # Optionally provide options files on the command line. Later files overwrite newer files
    for my $file (@ARGV){
        print "Reading options from $file\n";
        %$options = (%$options, %{Options::Read($file)});
    }
    
    AcceptableMCMC();
}

sub AcceptableMCMC {
    my $tree_filename = $options->{tree_filename};
    my $sequences_filename = $options->{sequences_filename};
    my $out_directory = $options->{out_directory};
    
    #$out_directory = UniqueDirectoryAppendNumber($out_directory);
    mkdir $out_directory if not -d $out_directory;
    Options::Write($options, $options->{out_directory} . "AcceptableMCMC_WAG.options");

    $options->{tree_out_filename} = $out_directory . $options->{tree_out_filename};
    $options->{sequences_out_filename} = $out_directory . $options->{sequences_out_filename};

    
    # MCMC files
    my $mcmc_out_filename = $out_directory . "mcmc";
    my $model_out_filename = $out_directory . "model";
    
    # data = [distance, C/D]
    print "Reading Data\n";
    my $data = ReadData($tree_filename, $sequences_filename);
    

    # Initialize output files
    my $mcmc_out = InitializeMcmcOutput($mcmc_out_filename);

    my $model = InitializeModel($options->{model_filename});


    print "Calculating WAG likelihood\n\n";
    

    my $aggregate = {
        node_sites => $number_of_branches * $number_of_sites,
        substitutions =>{},
        };
        
    #Tree::recurse_pre($data, \&Aggregator, $aggregate, $model);
    
    my $log_likelihood = LogLikelihoodOfModelGivenData($model, $data);
    RecordMcmcState($mcmc_out, 1, $log_likelihood, 1, 1, 0);
    
    print "Likelihood Calculation Finished\n";
}

sub LogLikelihoodNode {
    my ($node, $model, $log_likelihood) = @_;
    
    # sites are 1 indexed
    for my $site (1..$number_of_sites){
        LogLikelihoodNodeSite($node, $site,$model, $log_likelihood);
    }
}

sub LogLikelihoodNodeSite {
    my ($node, $site,$model, $log_likelihood) = @_;
    
    my $node_residue = substr($node->{Sequence}, $site-1, 1);
    
    if (not $node_residue){die "Node $node->{Name} has no sequence\n";}
    
    if ($node->{Up}) {
        my $anc_residue = substr($node->{Up}->{Sequence}, $site-1, 1);
        
        $$log_likelihood += log($model->{$anc_residue}->{$node_residue});
    }
    
    die "Log L is > 0 $node->{Name} $site: $$log_likelihood" if $$log_likelihood > 0;
}

    
sub LogLikelihoodOfModelGivenData {
    my ($model, $data) = @_;

    my $log_likelihood = 0;
    Tree::recurse_pre($data, \&LogLikelihoodNode, $model, \$log_likelihood);

    return $log_likelihood;
}
    

sub InitializeModel {
    my ($model_filename) = @_;
    
    open my $fh, "<", $model_filename or die "Can't find model file"; 
    
    my $exch = {};
    my $freqs = [];
    my $row_aa = 0;
    
    while (<$fh>){
        #print"line:$_\n";
        next if not $_;
        chomp;
        my $fields = [split " ", $_];
        if ($row_aa == 20){
            #print "@$fields";
            $freqs = $fields;
            last;
        }
        my $nfields = scalar @$fields;
         while (my ($index, $ex) = each(@$fields)){
             # These indices are wrong
             my $col_aa = $amino_acids->[$index+1];
             $exch->{$amino_acids->[$row_aa+1]}->{$amino_acids->[$index]} = $ex;
             $exch->{$amino_acids->[$index]}->{$amino_acids->[$row_aa+1]} = $ex;
         }
         $row_aa++;
    }
    print scalar@$freqs . ": @$freqs \n";
    for my $row_aa (@$amino_acids){
        while (my ($index, $col_aa) = each(@$amino_acids)){
            if ($row_aa ne $col_aa){
                # must multiply by max branch length to get the substitution probability
                $exch->{$row_aa}->{$col_aa} = $exch->{$row_aa}->{$col_aa} * $freqs->[$index] * $options->{max_branch_length};
            }
        }
    }
    for my $row_aa (@$amino_acids){
        $exch->{$row_aa}->{$row_aa} = 1;
        while (my ($index, $col_aa) = each(@$amino_acids)){
            if ($row_aa ne $col_aa){
                $exch->{$row_aa}->{$row_aa} -= $exch->{$row_aa}->{$col_aa};
            }
        }
    }
    
    my $model = $exch; 
    
    return $model;
}


sub ReadData {
    my ($tree_filename,$sequences_filename) = @_;
    
    my $tree = Tree::FromFile($tree_filename);
    
    my $sequences = Sequences::FromFasta($sequences_filename);
    Tree::AttachSequences($tree, $sequences);
    
    # a global
    $number_of_sites = length($tree->{Sequence});
    die "Number of sites could not be determined: $tree->{Name}" if not $number_of_sites;
    
    # B1 allows us to ignore double flip probability
    Tree::recurse_pre($tree, \&Tree::ToB1, $options->{max_branch_length});
    
    # In case the root has been segmented
    while (defined $tree->{Up}){ $tree = $tree->{Up} }
    
    $number_of_branches = scalar(@{Tree::GetAllNames($tree)});
    
    Tree::ToFile($tree, $options->{tree_out_filename});
    
    
    Tree::recurse_post($tree, \&PropagateSequences);
    
    my $new_sequences = Tree::ToSequences($tree);
    Sequences::ToFasta($new_sequences,$options->{sequences_out_filename});
    
    # globals for now
    $ancestral_node_names = Tree::GetAncestorNames($tree);
    $all_node_names = Tree::GetAllNames($tree);

    
    my $data = $tree;

    

    return $data;
}



sub PropagateSequences {
    my ($tree) = @_;
    
    if (not defined $tree->{Sequence}) {
        $tree->{Sequence} = $tree->{Left}->{Sequence};
    }
}


