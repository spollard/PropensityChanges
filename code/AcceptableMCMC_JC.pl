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
require "Progress.pl";
#require "UniqueDirectory.pl";


require "amino_acids.pl";
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
    $options = Options::Read("AcceptableMCMC_JC.options");

    # Optionally provide options files on the command line. Later files overwrite newer files
    for my $file (@ARGV){
        print "Reading options from $file\n";
        %$options = (%$options, %{Options::Read($file)});
    }
    
    AcceptableMCMC() if $options->{run_mcmc};
    #system("Rscript MakeMCMCGraphs.R $options->{out_directory} $options->{burnin}") if $options->{run_graphs};
    
}

sub AcceptableMCMC {
    my $tree_filename = $options->{tree_filename};
    my $sequences_filename = $options->{sequences_filename};
    my $generations = $options->{generations};
    my $out_directory = $options->{out_directory};
    
    #$out_directory = UniqueDirectoryAppendNumber($out_directory);
    mkdir $out_directory if not -d $out_directory;
    Options::Write($options, $options->{out_directory} . "AcceptableMCMC_JC.options");

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
    my $model_out = InitializeModelOutput($model_out_filename);

    my $model = InitializeModel();


    print "Beginning MCMC loop\n\n";
    
    my $progress = Progress($generations, $options->{time_between_prints});

    $empty_aggregate = {
        node_sites => $number_of_branches * $number_of_sites,
        substitutions => 0,
        };
        
    my $aggregate = clone($empty_aggregate);
    #Tree::recurse_pre($data, \&Aggregator, $aggregate, $model);
    

    # The next two parameters are rates and so the time needs to be accounted for
    
    #$aggregate->{optimal_sub_prob} = $aggregate->{substitutions} / $aggregate->{node_sites};
    
    my $proposed_log_likelihood = 0;
    my $log_likelihood = 0;
    if ($options->{use_aggregate}){
        $log_likelihood = LogLikelihoodOfModelGivenAggregate($model, $aggregate);
    } else {
        $log_likelihood = LogLikelihoodOfModelGivenData($model, $data);
    }

    open my $optimals, ">", $options->{out_directory} . "optimals";
    print $optimals join("\t", "optimal_acceptable_prior", "optimal_switch_prob", "optimal_sub_prob"),"\n";
    
    for my $generation (1 .. $generations) {
        $progress->($generation);
        
        if (not $options->{debug_gibbs}){
            my $proposed_model = clone($model);
            SampleProbabilities($proposed_model, $data);
            
            my $proposed_log_likelihood = 0;
            if ($options->{use_aggregate}){
                $proposed_log_likelihood = LogLikelihoodOfModelGivenAggregate($proposed_model, $aggregate);
                 if ($log_likelihood > 0){print Dumper $model, $aggregate; die "Aggregate log likelihood > 0";}
            } else {
                $proposed_log_likelihood = LogLikelihoodOfModelGivenData($proposed_model, $data);
                die "Non Aggregate log likelihood > 0" if $log_likelihood > 0;
            }
            
            if (log(rand) < $proposed_log_likelihood - $log_likelihood){
                $model = $proposed_model;
                $log_likelihood = $proposed_log_likelihood;
            }
        }
        
        RecordMcmcState($mcmc_out, $generation, $log_likelihood, 1, 1, 0);
        RecordModelState($model, $model_out);
    }
    print "MCMC Finished\n";
}


sub InitializeModelOutput {
    my ($model_out_filename, $amino_acids_considered, $hidden_states_out_filename) = @_;
    open my $model_out, ">", $model_out_filename;
    
    print $model_out "SubRate\n";
    
    return $model_out;
}

sub RecordModelState {
    my ($model, $model_out) = @_;
    
    print $model_out $model->{SubRate}, "\n";
    

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
        
        my $node_SubRate = 1;
        # residue sub
        if ($anc_residue ne $node_residue) {
            # Substitution
            $node_SubRate = $model->{SubRate}#*$node->{Distance};
         } else {
             # No substitution
             $node_SubRate = 1-$model->{SubRate}#*$node->{Distance};
         }
        
        
        $$log_likelihood += log($node_SubRate);
      
    }
    die "Log L is > 0 $node->{Name} $site: $$log_likelihood" if $$log_likelihood > 0;
}


    
sub LogLikelihoodOfModelGivenData {
    my ($model, $data) = @_;

    my $log_likelihood = 0;
    Tree::recurse_pre($data, \&LogLikelihoodNode, $model, \$log_likelihood);

    return $log_likelihood;
}
    

    
sub SampleProbabilities {
    my ($model, $aggregate) = @_;

    $model->{SubRate} = Samplers::normal_p($model->{SubRate}, $options->{sub_stdev}) if $options->{sample_sub_rate};
}


sub InitializeModel {
    my () = @_;
    my $model = {
        SubRate => $options->{SubRate} // rand,
    };
    
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
    die "Number of branches could not be determined" if not $number_of_branches;
    
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
