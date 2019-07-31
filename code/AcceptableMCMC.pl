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
    $options = Options::Read("AcceptableMCMC.options");

    # Optionally provide options files on the command line. Later files overwrite newer files
    for my $file (@ARGV){
        print "Reading options from $file\n";
        %$options = (%$options, %{Options::Read($file)});
    }
    
    AcceptableMCMC() if $options->{run_mcmc};
    system("Rscript MakeMCMCGraphs.R $options->{out_directory} $options->{burnin}") if $options->{run_graphs};
    system("perl stats.pl @ARGV") if $options->{run_stats};
    system("perl collect_probable_acceptabilities.pl $options->{out_directory}/sites $options->{out_directory}/acceptabilities2.fasta") if $options->{run_stats};
}

sub AcceptableMCMC {
    my $tree_filename = $options->{tree_filename};
    my $sequences_filename = $options->{sequences_filename};
    my $generations = $options->{generations};
    my $out_directory = $options->{out_directory};
    
    #$out_directory = UniqueDirectoryAppendNumber($out_directory);
    mkdir $out_directory if not -d $out_directory;
    Options::Write($options, $options->{out_directory} . "AcceptableMCMC.options");

    $options->{tree_out_filename} = $out_directory . $options->{tree_out_filename};
    $options->{sequences_out_filename} = $out_directory . $options->{sequences_out_filename};
    $options->{initial_acceptables_out_filename} = $out_directory . $options->{initial_acceptables_out_filename};

    $options->{sites_dir} = $out_directory . "sites";
    mkdir $options->{sites_dir} if not -d $options->{sites_dir};

    # MCMC files
    my $mcmc_out_filename = $out_directory . "mcmc";
    my $model_out_filename = $out_directory . "model";
    my $hidden_states_out_filename = $out_directory . "hidden_states";

    # data = [distance, C/D]
    print "Reading Data\n";
    my $data = ReadData($tree_filename, $sequences_filename);
    

    # Initialize output files
    my $mcmc_out = InitializeMcmcOutput($mcmc_out_filename);
    my ($model_out, $hidden_states_out) = InitializeModelOutput($model_out_filename, $amino_acids_considered,$hidden_states_out_filename);

    my $model = InitializeModel($amino_acids_considered, $data);


    print "Beginning MCMC loop\n\n";
    
    my $progress = Progress($generations, $options->{time_between_prints});

    $empty_aggregate = {
        # subtract 1 bc the root cant sub or switch
        possible_switches => ($number_of_branches-1) * $number_of_amino_acids * $number_of_sites,
        node_sites => ($number_of_branches-1) * $number_of_sites,
        number_of_substitutions => 0,
        substitutions => {},
        transitions => 0,
        transversions => 0,
        number_of_switches => 0,
        subs_q_a => 0,
        subs_q_u => 0,
        sites_with_acceptable_residue =>0,
        acceptables_sizes =>[],
        };
        
    my $aggregate = clone($empty_aggregate);
    Tree::recurse_pre($data, \&Aggregator, $aggregate, $model);
    $aggregate->{optimal_acceptable_prior} = $aggregate->{sites_with_acceptable_residue} / $number_of_branches / $number_of_sites;

    # The next two parameters are rates and so the time needs to be accounted for
    $aggregate->{optimal_switch_prob} = $aggregate->{number_of_switches} / $aggregate->{possible_switches};
    $aggregate->{optimal_sub_prob} = $aggregate->{number_of_substitutions} / $aggregate->{node_sites};
    
    if ($options->{save_aggregate}) {
        SaveAggregate($aggregate);exit;
    }
    
    my $proposed_log_likelihood = 0;
    my $log_likelihood = 0;;
    if ($options->{use_aggregate}){
        $log_likelihood = LogLikelihoodOfModelGivenAggregate($model, $aggregate);
    } else {
        $log_likelihood = LogLikelihoodOfModelGivenData($model, $data);
    }

    open my $optimals, ">", $options->{out_directory} . "optimals";
    print $optimals join("\t", "optimal_acceptable_prior", "optimal_switch_prob", "optimal_sub_prob"),"\n";
    
    for my $generation (1 .. $generations) {
        $progress->($generation);
        if ($options->{sample_acceptabilities} and rand() < $options->{resample_acceptabilities_probability}) {
            if ($options->{gibbs_aggregate}) {
                Tree::recurse_pre($data, \&GibbsSampleAcceptablesViaAggregation, $model);
            } else {
                Tree::recurse_pre($data, \&GibbsSampleAcceptables, $model);
            }
            if ($options->{use_aggregate}){
                # Reset the aggregator
                $aggregate = clone($empty_aggregate);
                Tree::recurse_pre($data, \&Aggregator, $aggregate, $model);
                #print(Dumper $aggregate);<STDIN>;
                # Calculate optima
                $aggregate->{optimal_acceptable_prior} = $aggregate->{sites_with_acceptable_residue} / $number_of_branches / $number_of_sites;
        
                # The next two parameters are rates and so the time needs to be accounted for
                $aggregate->{optimal_switch_prob} = $aggregate->{number_of_switches} / $aggregate->{possible_switches};
                $aggregate->{optimal_sub_prob} = $aggregate->{number_of_substitutions} / $aggregate->{node_sites};
            }
            # Could print optimals here, which is only when they change and not every generation
            #print $optimals  join("\t", $aggregate->{optimal_acceptable_prior}, $aggregate->{optimal_switch_prob}, $aggregate->{optimal_sub_prob}),"\n";
        }
        
        print $optimals  join("\t", $aggregate->{optimal_acceptable_prior}, $aggregate->{optimal_switch_prob}, $aggregate->{optimal_sub_prob}),"\n";
        
        if (not $options->{debug_gibbs}){
            my $proposed_model = clone($model);
            SampleProbabilities($proposed_model, $data);
            
            my $proposed_log_likelihood = 0;
            if ($options->{use_aggregate}){
                die "Cannot use aggregate when using jukes-cantor or wag" if $options->{use_jc} or $options->{use_wag};
                $proposed_log_likelihood = LogLikelihoodOfModelGivenAggregate($proposed_model, $aggregate);
                #$log_likelihood = LogLikelihoodOfModelGivenAggregate($model, $aggregate);
                 if ($log_likelihood > 0){print Dumper $model, $aggregate; die "Aggregate log likelihood > 0";}
            } else {
                $proposed_log_likelihood = LogLikelihoodOfModelGivenData($proposed_model, $data);
                #$log_likelihood = LogLikelihoodOfModelGivenData($model, $data);
                die "Non Aggregate log likelihood > 0" if $log_likelihood > 0;
            }
            
    #        print STDERR join("\t",$log_likelihood_d, $log_likelihood), "\n";
            
            if (log(rand) < $proposed_log_likelihood - $log_likelihood){
                $model = $proposed_model;
                $log_likelihood = $proposed_log_likelihood;
            }
        }
        
        RecordMcmcState($mcmc_out, $generation, $log_likelihood, 1, 1, 0);
        RecordModelState($model, $model_out);
        
        Tree::recurse_pre($data, \&RecordAcceptables, $hidden_states_out);
        print $hidden_states_out "\n";
        
        RecordAcceptablesToSites($data);
    }
    print "MCMC Finished\n";
}


sub InitializeModelOutput {
    my ($model_out_filename, $amino_acids_considered, $hidden_states_out_filename) = @_;
    open my $model_out, ">", $model_out_filename;
    
    print $model_out join("\t", qw(Prior_A SwitchRate SubRate FewAcceptables)), "\n";

    for my $site (1..$number_of_sites){
        open my $fh2, ">", $options->{sites_dir} . "/" . $site;
        print $fh2 join("\t", @$all_node_names),"\n";
    }
    
    open my $hidden_states_out, ">", $hidden_states_out_filename;
    
    return $model_out,$hidden_states_out;
}

sub RecordModelState {
    my ($model, $model_out) = @_;
    
    print $model_out join("\t", ($model->{Prior_A}, $model->{SwitchProb}, $model->{SubProb}, $model->{FewAcceptables})), "\n";
    

}
sub RecordAcceptables   {
    my ($node, $hidden_states_out) = @_;
    print $hidden_states_out ">$node->{Name}\n" . join("|", map({join"",$node->{Acceptables}->[$_]->members()} (0..($number_of_sites-1))) ) . "\n";
}

sub RecordAcceptablesToSites   {
    my ($data) = @_;
        my $acceptables = {};
        Tree::recurse_pre($data, \&CollectAcceptables2, $acceptables);
        
        for my $site (1..$number_of_sites){
            open my $fh2, ">>", $options->{sites_dir} . "/" . $site;
            print $fh2 join("\t", map({$acceptables->{$_}->[$site-1]->as_string()} @$all_node_names) ) . "\n";
        }
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
    
    #my $node_acceptables = $model->{Acceptables}->{$node->{Name}}->[$site-1];
    my $node_acceptables = $node->{Acceptables}->[$site-1];
    my $node_residue_acceptable = $node_acceptables->has($node_residue)? 1 : 0;
    
    
    # Prior of being in acceptable or not
    # May be unnecessary if the probability of substituting to an unacceptable 
    # aa is low enough
    $$log_likelihood += log($node_acceptables->has($node_residue) ? $model->{Prior_A} : 1-$model->{Prior_A});
    
    # Prior against having a high number of acceptables
    # Might not be necessary if the sub rate depends on the number 
    # of alternative possible substitutions
    $$log_likelihood += log(ProbOfNumberAcceptables($model,$node_acceptables->size()));
    
    if ($node->{Up}) {
        my $anc_residue = substr($node->{Up}->{Sequence}, $site-1, 1);
        my $anc_acceptables = $node->{Up}->{Acceptables}->[$site-1];
        my $anc_residue_acceptable = $anc_acceptables->has($anc_residue)? 1 : 0;
        
        my $node_subprob = 1;
        # residue sub
        if ($anc_residue ne $node_residue) {
            # Substitution
            if ($anc_residue_acceptable and $node_residue_acceptable){
                # Need transcriptions and transversions here
                # number of sub options
                my $sub_options = $node_acceptables->size() > 1 ? $node_acceptables->size()-1 : 1;
                # This will help keep the acceptable sets small
                $node_subprob = $model->{SubProb} / $sub_options;
                
                if (transition($anc_residue, $node_residue)){
                    #$node_subprob = $model->{TransitionProb} / $sub_options;
                } else{
                    #$node_subprob = $model->{TransversionProb} / $sub_options;
                }
            } elsif (not $anc_residue_acceptable and $node_residue_acceptable){
                 $node_subprob = $model->{q_a};
             } elsif ($anc_residue_acceptable and not $node_residue_acceptable){
                 # If q_u is very low, then the resident aa will be 
                 # acceptale often
                 $node_subprob = $model->{q_u};
             } else { # both not acceptable
                 # If q_u_u is very low, then the resident aa will be 
                 # acceptale often
                $node_subprob = $model->{q_u_u};
             }
         } else {
             # No substitution
             $node_subprob = 1-$model->{SubProb};
         }
        
        die "Prob is negative: $node_subprob; subprob $model->{SubProb}; acc size: ".$node_acceptables->size() if $node_subprob < 0;
        
        # Since we don't know if the sub happened on the top half or bottom half of the branch
        #$$log_likelihood += log(0.5 * $anc_subprob + 0.5 * $node_subprob);
        $$log_likelihood += log($node_subprob);
        
        # Prior of anc being in acceptable or not
        #$$log_likelihood += log($node_acceptables->has($anc_residue) ? $model->{Prior_A} : 1-$model->{Prior_A});
        
        # How many residues switched to/from acceptable
        my $differences = $anc_acceptables->symmetric_difference($node_acceptables)->size();
        my $sames = $number_of_amino_acids - $differences;
        
        if ($differences){
            $$log_likelihood += $differences*log($model->{SwitchProb});
        }
        if ($sames){
            $$log_likelihood += $sames * log(1-$model->{SwitchProb});
        }
    }
    
    die "Log L is > 0 $node->{Name} $site: $$log_likelihood" if $$log_likelihood > 0;
}

sub ProbOfNumberAcceptables {
    my ($model,$number_of_acceptables) = @_;
    #return $model->{FewAcceptables}**$number_of_acceptables;
    
    return 0.00001 if $number_of_acceptables > $#{$options->{priors}};
    
    return $options->{priors}->[$number_of_acceptables];
}
    
sub LogLikelihoodOfModelGivenData {
    my ($model, $data) = @_;

    my $log_likelihood = 0;
    Tree::recurse_pre($data, \&LogLikelihoodNode, $model, \$log_likelihood);

    return $log_likelihood;
}
    
sub LogLikelihoodOfModelGivenAggregate {
    my ($model, $aggregate) = @_;
    
        
    if ($options->{debug_aggregate_likelihood}){
        print Dumper $aggregate;
        print "Log L from switches: " . $aggregate->{number_of_switches} * log($model->{SwitchProb}), "\n";
        print "Log L from not switches: " . ($aggregate->{possible_switches} - $aggregate->{number_of_switches}) * log(1-$model->{SwitchProb}), "\n";
        print "Log L from subs: " . $aggregate->{number_of_substitutions} * log($model->{SubProb}), "\n";
        print "Log L from not subs: " . ($aggregate->{node_sites} - $aggregate->{number_of_substitutions}) * log(1-$model->{SubProb}), "\n";
        print "Log L from acceptable residues: " . $aggregate->{sites_with_acceptable_residue} * log($model->{Prior_A}), "\n";
        print "Log L from not acceptable residues: " . ($aggregate->{node_sites} - $aggregate->{sites_with_acceptable_residue}) * log(1-$model->{Prior_A}), "\n";
        print "priors: @{$options->{priors}}\n";
        for my $i (0..$#{$aggregate->{acceptables_sizes}}){
            print "Log L from $i acceptables: " . (($aggregate->{acceptables_sizes}->[$i]//0) * log($options->{priors}->[$i])), "\n";
        }
    }
    
    my $log_likelihood = 0;
    #log likelihood from switches and not switches
    my $switch_likelihood = $aggregate->{number_of_switches} * log($model->{SwitchProb});
    die "Switch likelihood > 0" if $switch_likelihood > 0;
    my $not_switch_likelihood = ($aggregate->{possible_switches} - $aggregate->{number_of_switches}) * log(1-$model->{SwitchProb});
    die "Not switch likelihood > 0" if $not_switch_likelihood > 0;
    my $sub_likelihood = $aggregate->{number_of_substitutions} * log($model->{SubProb});
    die "Sub likelihood > 0" if $sub_likelihood > 0;
    my $not_sub_likelihood = ($aggregate->{node_sites} - $aggregate->{number_of_substitutions}) * log(1-$model->{SubProb});
    die "Not sub likelihood > 0" if $not_sub_likelihood > 0;
    my $acceptable_res_likelihood = $aggregate->{sites_with_acceptable_residue} * log($model->{Prior_A});
    die "Acceptable res likelihood > 0" if $acceptable_res_likelihood > 0;
    my $unacceptable_res_likelihood = ($aggregate->{node_sites} - $aggregate->{sites_with_acceptable_residue}) * log(1-$model->{Prior_A});
    die "Unacceptable res likelihood > 0" if $unacceptable_res_likelihood > 0;
    
    my $sizes_likelihoods = [];
    for my $i (0..$#{$aggregate->{acceptables_sizes}}){
        my $size_like = (($aggregate->{acceptables_sizes}->[$i]//0) * log(ProbOfNumberAcceptables($model, $i)));
        die "Size likelihood for size $i > 0" if $size_like > 0;
        push @$sizes_likelihoods,$size_like;
    }
    
    $log_likelihood  = $switch_likelihood + $not_switch_likelihood + $sub_likelihood + $not_sub_likelihood + $acceptable_res_likelihood + $unacceptable_res_likelihood;
    for (@$sizes_likelihoods){
        $log_likelihood += $_;
    }
    #~ $log_likelihood += $aggregate->{number_of_switches} * log($model->{SwitchProb}) + ($aggregate->{possible_switches} - $aggregate->{number_of_switches}) * log(1-$model->{SwitchProb});
    
    #~ # log likelihood from subs and not subs
    #~ $log_likelihood += $aggregate->{number_of_substitutions} * log($model->{SubProb}) + ($aggregate->{node_sites} - $aggregate->{number_of_substitutions}) * log(1-$model->{SubProb});
    
    #~ $log_likelihood += $aggregate->{sites_with_acceptable_residue} * log($model->{Prior_A}) + ($aggregate->{node_sites} - $aggregate->{sites_with_acceptable_residue}) * log(1-$model->{Prior_A});
    
    #~ for my $i (0..$#{$aggregate->{acceptables_sizes}}){
        #~ $log_likelihood += ($aggregate->{acceptables_sizes}->[$i]//0) * log(ProbOfNumberAcceptables($model, $i));
    #~ }

    
    return $log_likelihood
}


sub GibbsSampleAcceptables {
    my ($node, $model) = @_;
    print "Gibbs sampling acceptables" if $options->{debug_gibbs};
    for my $site (1..$number_of_sites){
        #my $acceptables = $model->{Acceptables}->{$node->{Name}}->[$site-1];
        my $acceptables = $node->{Acceptables}->[$site-1];
        # Need to loop over residues in a random order
        for my $residue (shuffle @$amino_acids_considered) {
            print "Testing if $residue is acceptable at $site\n" if $options->{debug_gibbs};
            
            my $aggregate = clone($empty_aggregate);
            $aggregate->{node_sites} = 1;
            $aggregate->{possible_switches} = 1 * $number_of_amino_acids;
            
            $acceptables->insert($residue);
            my $log_likelihood_acceptable = 0;
            LogLikelihoodNodeSite($node, $site, $model, \$log_likelihood_acceptable);
            LogLikelihoodNodeSite($node->{Left}, $site, $model, \$log_likelihood_acceptable) if $node->{Left};
            LogLikelihoodNodeSite($node->{Right}, $site, $model, \$log_likelihood_acceptable) if $node->{Right};
            
            
            $acceptables->remove($residue);
            
            my $log_likelihood_unacceptable = 0;
            LogLikelihoodNodeSite($node, $site, $model, \$log_likelihood_unacceptable);
            LogLikelihoodNodeSite($node->{Left}, $site, $model, \$log_likelihood_unacceptable) if $node->{Left};
            LogLikelihoodNodeSite($node->{Right}, $site, $model, \$log_likelihood_unacceptable) if $node->{Right};
            
            #print "prob $residue Acceptable $log_likelihood_acceptable; prob not $log_likelihood_unacceptable\n"  . exp($log_likelihood_acceptable) . "      " . exp($log_likelihood_unacceptable) . "\n" if $options->{debug_gibbs};
            
            # These are log probabilities
            print "Prob of $residue being acceptable at $node->{Name} $site:" . exp($log_likelihood_acceptable)/(exp($log_likelihood_unacceptable)+exp($log_likelihood_acceptable)) . "\n\n" if $options->{debug_gibbs};
            if (rand() < exp($log_likelihood_acceptable)/(exp($log_likelihood_unacceptable)+exp($log_likelihood_acceptable))) {
                $acceptables->insert($residue);
            }
            #else {
                # It has already been removed so you don't have to do this
                #$acceptables->remove($residue);
            #}
            
        }
    }
}


sub GibbsSampleAcceptablesViaAggregation {
    my ($node, $model) = @_;
    print "Gibbs sampling acceptables" if $options->{debug_gibbs};
    for my $site (1..$number_of_sites){
        my $acceptables = $node->{Acceptables}->[$site-1];
        # Need to loop over residues in a random order
        for my $residue (shuffle @$amino_acids_considered) {
            print "Testing if $residue is acceptable at $site\n" if $options->{debug_gibbs};
            
            my $aggregate = clone($empty_aggregate);
            $aggregate->{node_sites} = 1;
            $aggregate->{possible_switches} = 1 * $number_of_amino_acids;
            
            $acceptables->insert($residue);
            
            AggregatorNodeSite($node, $aggregate, $model, $site);
            if ($node->{Left}){
                AggregatorNodeSite($node->{Left}, $aggregate, $model, $site);
                $aggregate->{node_sites} += 1;
                $aggregate->{possible_switches} += 1 * $number_of_amino_acids;
            }
            if ($node->{Right}){
                AggregatorNodeSite($node->{Right}, $aggregate, $model, $site);
                $aggregate->{node_sites} += 1;
                $aggregate->{possible_switches} += 1 * $number_of_amino_acids;
            }
            
            my $log_likelihood_acceptable = LogLikelihoodOfModelGivenAggregate($model, $aggregate);
            print "Aggregate if $residue is acceptable\n" if $options->{debug_gibbs};
            print Dumper $aggregate if $options->{debug_gibbs};
            
            
            # reset the aggregate
            $aggregate = clone($empty_aggregate);
            $aggregate->{node_sites} = 1;
            $aggregate->{possible_switches} = 1 * $number_of_amino_acids;
            
            
            $acceptables->remove($residue);
            
            AggregatorNodeSite($node, $aggregate, $model, $site);
            if ($node->{Left}){
                AggregatorNodeSite($node->{Left}, $aggregate, $model, $site);
                $aggregate->{node_sites} += 1;
                $aggregate->{possible_switches} += 1 * $number_of_amino_acids;
                }
            if ($node->{Right}){
                AggregatorNodeSite($node->{Right}, $aggregate, $model, $site);
                $aggregate->{node_sites} += 1;
                $aggregate->{possible_switches} += 1 * $number_of_amino_acids;
                }
            
            print "Aggregate if $residue is unacceptable\n" if $options->{debug_gibbs};
            print Dumper $aggregate if $options->{debug_gibbs};
            
            my $log_likelihood_unacceptable = LogLikelihoodOfModelGivenAggregate($model, $aggregate);
            
            #print "prob $residue Acceptable $log_likelihood_acceptable; prob not $log_likelihood_unacceptable\n"  . exp($log_likelihood_acceptable) . "      " . exp($log_likelihood_unacceptable) . "\n" if $options->{debug_gibbs};
            
            # These are log probabilities
            print "Prob of $residue being acceptable at $node->{Name} $site:" . exp($log_likelihood_acceptable)/(exp($log_likelihood_unacceptable)+exp($log_likelihood_acceptable)) . "\n\n" if $options->{debug_gibbs};
            if (rand() < exp($log_likelihood_acceptable)/(exp($log_likelihood_unacceptable)+exp($log_likelihood_acceptable))) {
                $acceptables->insert($residue);
            }
            #else {
                # It has already been removed so you don't have to do this
                #$acceptables->remove($residue);
            #}
            
        }
    }
}
    
sub SampleProbabilities {
    my ($model, $aggregate) = @_;

    $model->{Prior_A} = Samplers::normal_p($model->{Prior_A},$options->{prior_stdev}) if $options->{sample_residue_acceptability_prior};
    if ($model->{Prior_A} < 0.5){
        $model->{Prior_A} = 0.6;
    }
    
    $model->{SwitchProb} = Samplers::normal_p($model->{SwitchProb}, $options->{switch_stdev}) if $options->{sample_switch_rate};
    $model->{SubProb} = Samplers::normal_p($model->{SubProb}, $options->{sub_stdev}) if $options->{sample_sub_rate};
}


sub Aggregator {
    my ($node, $aggregate, $model) = @_;
    for my $site (1..$number_of_sites){
        AggregatorNodeSite($node, $aggregate, $model, $site);
    }
}
        
sub AggregatorNodeSite {
    my ($node, $aggregate, $model, $site) = @_; 

    my $node_acceptables = $node->{Acceptables}->[$site-1];
    my $node_residue = substr($node->{Sequence}, $site-1, 1);
    
    my $node_residue_acceptable = $node_acceptables->has($node_residue) ? 1 : 0;
    
    $aggregate->{sites_with_acceptable_residue}++ if $node_residue_acceptable;
    $aggregate->{acceptables_sizes}->[$node_acceptables->size()]++;
    
    if ($node->{Up}) {
        my $anc_residue = substr($node->{Up}->{Sequence}, $site-1, 1);
        my $anc_acceptables = $node->{Up}->{Acceptables}->[$site-1];
        
        my $anc_residue_acceptable = $anc_acceptables->has($anc_residue)? 1 : 0;
        
        $aggregate->{substitutions}->{$anc_residue}->{$node_residue}++;
        
         if ($anc_residue ne $node_residue){
             $aggregate->{number_of_substitutions}++;
            # Need transcriptions and transversions here
            if ($anc_residue_acceptable and $node_residue_acceptable){
                
                # This does notwork for amino acids yet
                if (transition($anc_residue, $node_residue)){
                    $aggregate->{transitions}++;
                } else{
                    $aggregate->{transversions}++;
                }
            } elsif (not $anc_residue_acceptable and $node_residue_acceptable){
             $aggregate->{q_a}++;
            } elsif ($anc_residue_acceptable and not $node_residue_acceptable){
             $aggregate->{q_u}++;
            } else { #both not acceptable
             $aggregate->{q_u_u}++;
            }
         }
        
        $aggregate->{number_of_switches} += $anc_acceptables->symmetric_difference($node_acceptables)->size();
    }
}

sub SaveAggregate {
    my ($aggregate) = @_;
    
    open my $fh, ">", $options->{out_directory} . "aggregate" or die;
    print $fh Dumper $aggregate;exit;
    print $fh "Number of substitutions: $aggregate->{number_of_substitutions}\n";
    print $fh "Number of switches: $aggregate->{number_of_switches}\n";
    my $string = join("\t", @$amino_acids). "\n";
    for my $row_aa (@$amino_acids){
        while (my ($index, $col_aa) = each(@$amino_acids)){
                $string .= ($index == 0 ? "" : "\t") . ($aggregate->{substitutions}->{$row_aa}->{$col_aa}//0);
        }
        $string .= "\n";
    }
    
}


sub transition {
    my ($aa1, $aa2) = @_;
    
    return $aa2 eq "T" if $aa1 eq "A";
    return $aa2 eq "G" if $aa1 eq "C";
}


sub InitializeModel {
    my ($amino_acids_considered, $data) = @_;
    my $model = {
        # The following parameters could be sampled, but I'd rather 
        # have them be input by the user. This also reduces the 
        # amount of sampling required.
        # Priors for a site to have a residue which is Acceptable
        Prior_A => $options->{Prior_A} // (1 - 0.5 * rand), # must be > 0.5
        # why penalize substitutions?
        SubProb => $options->{SubProb} // rand,
        # Prob of an amino acid switching to/from acceptable
        SwitchProb => $options->{SwitchProb} // rand, 
        # Prior for a residue to be acceptable
        # does nothing at the moment
        FewAcceptables => rand,
        q_a => $options->{q_a} // rand,
        q_u => $options->{q_u} // rand,
        q_u_u => $options->{q_u_u} // rand,
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
    
    if ($options->{acceptables_filename}){
        my $acceptables = Sequences::FromFasta($options->{acceptables_filename});
        while (my ($name, $acc) = each(%$acceptables)){
            $acceptables->{$name} = [split /\|/, $acc, -1];
            die "Not enough acceptable sets provided" if not scalar @{$acceptables->{$name}} == $number_of_sites;
            for my $site (1..$number_of_sites){
                $acceptables->{$name}->[$site-1] = Set::Tiny->new(split "", $acceptables->{$name}->[$site-1]);
            }
        }
        Tree::recurse_pre($tree, \&InitializeAcceptables, $acceptables);
    } else {
        Tree::recurse_pre($tree, \&InitializeAcceptables);
    }
    
    # B1 allows us to ignore double flip probability
    Tree::recurse_pre($tree, \&Tree::ToB1, $options->{max_branch_length});
    
    # In case the root has been segmented
    while (defined $tree->{Up}){ $tree = $tree->{Up} }
    
    $number_of_branches = scalar(@{Tree::GetAllNames($tree)});
    
    Tree::ToFile($tree, $options->{tree_out_filename});
    
    
    Tree::recurse_post($tree, \&PropagateSequences);
    Tree::recurse_post($tree, \&PropagateAcceptables);
     
    my $new_sequences = Tree::ToSequences($tree);
    Sequences::ToFasta($new_sequences,$options->{sequences_out_filename});
    
    my $acceptables = {};
    Tree::recurse_pre($tree, \&CollectAcceptables, $acceptables);
    Sequences::ToFasta($acceptables,$options->{initial_acceptables_out_filename});

    # globals for now
    $ancestral_node_names = Tree::GetAncestorNames($tree);
    $all_node_names = Tree::GetAllNames($tree);

    
    # Collect numbers of events in the data
    
    
    
    my $data = $tree;

    

    return $data;
}



sub PropagateSequences {
    my ($tree) = @_;
    
    if (not defined $tree->{Sequence}) {
        $tree->{Sequence} = $tree->{Left}->{Sequence};
    }
}

sub PropagateAcceptables {
    my ($tree) = @_;
    
    if (not defined $tree->{Acceptables}) {
        $tree->{Acceptables} = clone($tree->{Left}->{Acceptables});
    }
}

sub CollectAcceptables {
    my ($tree, $acceptables) = @_;
    
    $acceptables->{$tree->{Name}} = join("|", map({join "",$tree->{Acceptables}->[$_]->members()} (0..($number_of_sites-1))) )
}

sub CollectAcceptables2 {
    my ($tree, $acceptables) = @_;
    
    $acceptables->{$tree->{Name}} = $tree->{Acceptables};
}

sub InitializeAcceptables {
    my ($tree, $acceptables) = @_;
    if ($acceptables and exists $acceptables->{$tree->{Name}}) {
        $tree->{Acceptables} = $acceptables->{$tree->{Name}};
    } else {
        for my $site (1..$number_of_sites){
            $tree->{Acceptables}->[$site-1] = Set::Tiny->new(Samplers::random_choice($amino_acids_considered));
        }
    }
}



sub PrintAcceptables {
    my ($model) = @_;
    for my $site (1..$number_of_sites){
        print join("\t", map({$_ . ": " . $model->{Acceptables}->{$_}->[$site-1]->as_string()} @$all_node_names) ) . "\n";
    }
}
