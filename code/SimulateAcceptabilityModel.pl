# Simulate an Acceptablity model over a tree
#
#  I THINK THIS IS INCOMPLETE AND SIMPLY A COPY OF THE SimulateSingleAcceptabilityModel.pl
# Stephen Pollard 2019-3-4

use Math::Random;
use Data::Dumper;

use List::Util 'first';
use Set::Tiny;

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;


use Tree;
use Samplers;
use Propensities;
use Options;

require "amino_acids.pl";
our $amino_acids;


# Every node in a tree must be in a state. Along a branch, the state can
# change, just like amino acids. This makes sense because a node is an
# instant in time and so there are no state changes there.

unless (caller) {
    
     # Default options
    my $options = Options::Read("SimulateSingleAcceptabilityModel.options");
    
    # Provide options files on the command line. Later files overwrite newer files
    for my $file (@ARGV){
        print("Reading options from $file\n");
        %$options = (%$options, %{Options::Read($file)});
    }
    
    SimulateSingleAcceptabilityModel($options);
}

sub SimulateSingleAcceptabilityModel {
    my ($options) = @_;
    
    print(Dumper($options));
    
    my $model = {
        switch_rate => $options->{switch_rate},
        sub_rate => $options->{sub_rate},
        q_a => $options->{q_a},
        q_u => $options->{q_u},
        residues => $amino_acids,
    };
    
    open my $fh, "<", $options->{tree_filename};
    my $string = scalar <$fh>;
    
    #$string = Tree::GenerateNamesForInternalNodes($string);
    
    my $tree = Tree::FromString($string);
    
    
    #~ my $long_branches = Tree::FindLongBranches($tree, $options->{max_branch_length});
    #~ if (scalar @$long_branches) {
        #~ die "Tree has branches which are too long. Please break up branches until all are shorter than $options->{max_branch_length}"
    #~ }
    Tree::recurse_pre($tree, \&Tree::ToB1, $options->{max_branch_length});
    Tree::ToFile($tree, $options->{tree_filename} . ".out");

    # Initial acceptable amino acid
    my $aa = Samplers::random_choice($amino_acids);
    $tree->{HiddenState} = Set::Tiny->new($aa);
    $tree->{Sequence} = $aa;

    Tree::recurse_pre($tree, \&SimulateAcceptability, $model);

    open my $sequences_file, ">", $options->{sequences_filename};
    print $sequences_file Tree::AllSequencesToString($tree);

    open my $hiddenStates_file, ">", $options->{hidden_states_filename};
    print $hiddenStates_file Tree::AllHiddenStatesToString2($tree);
}


sub SimulateAcceptability {
    my ($descendant, $model) = @_;

#    print "Working with tree " . $descendant->{Name}, "\n";
    
    # Root does not have an ancestor
    return if not exists $descendant->{Up};
    
    my $ancestor = $descendant->{Up};
    
    $descendant->{HiddenState} = $ancestor->{HiddenState}->copy();

# Sample descendant acceptabilities
    for my $aa (@{$model->{residues}}) {
        if (rand() < $model->{switch_rate}*$descendant->{Distance}) {
            if ($descendant->{HiddenState}->member($aa)){
                $descendant->{HiddenState}->remove($aa);
            }
            else{
                $descendant->{HiddenState}->insert($aa);
            }
        }
    }
    
    # Adjust the size to the expected distribution
    # Without this there is no penalty to adding more and more 
    # amino acids to the acceptable set
    my $set_size = Samplers::weighted_random_choice([0.985, 0.015, 0.0001])+1;
    while ($descendant->{HiddenState}->size() > $set_size){
        $descendant->{HiddenState}->remove(Samplers::random_choice([$descendant->{HiddenState}->members()]));
    }
    while ($descendant->{HiddenState}->size() < $set_size){
        $descendant->{HiddenState}->insert(Samplers::random_choice($model->{residues}));
    }
    
    if (not $descendant->{HiddenState}->member($ancestor->{Sequence})){
        #$descendant->{HiddenState}->remove(Samplers::random_choice([$descendant->{HiddenState}->members()]));
        $descendant->{HiddenState}->insert($ancestor->{Sequence});
    }
    
# Then sample the descendant state
    my $probs = [];
    my $anc_was_acceptable = $ancestor->{HiddenState}->member($ancestor->{Sequence}) ? 1 : 0;
    for my $aa (@{$model->{residues}}) {
        
        my $dec_is_acceptable = $descendant->{HiddenState}->member($aa) ? 1 : 0;
        # Prior for the amino acid to be acceptable
        my $prob = $dec_is_acceptable ? 0.99 : 0.01;
        if ($aa eq $ancestor->{Sequence}) {
            # no substitution
            if ($anc_was_acceptable and $dec_is_acceptable){
                $prob *= 1 - $model->{sub_rate} * $descendant->{Distance};
            } elsif (not $anc_was_acceptable and $dec_is_acceptable) {
                $prob *= 1 - $model->{q_a} * $descendant->{Distance};
            }elsif ($anc_was_acceptable and not $dec_is_acceptable) {
                $prob *= 1 - $model->{q_u} * $descendant->{Distance};
            } else { # Both not acceptable
                $prob *= 1 - $model->{sub_rate} * $descendant->{Distance};
            }
        }
        else {
            # Substitution
            if ($anc_was_acceptable and $dec_is_acceptable){
                $prob *= $model->{sub_rate} * $descendant->{Distance};
            } elsif (not $anc_was_acceptable and $dec_is_acceptable) {
                $prob *= $model->{q_a} * $descendant->{Distance};
            } elsif ($anc_was_acceptable and not $dec_is_acceptable) {
                $prob *= $model->{q_u} * $descendant->{Distance};
            } else { # Both not acceptable
                $prob *= $model->{sub_rate} * $descendant->{Distance};
            }
        }
        push(@$probs, $prob);
    }
    #~ print $ancestor->{HiddenState}->as_string(), "\t", $ancestor->{Sequence}, "\n";
    #~ print $descendant->{HiddenState}->as_string(), "\n";
    #~ for my $i (0 .. $#{$model->{residues}}){
        #~ print("$model->{residues}->[$i]: $probs->[$i]\t");
    #~ }
    #~ print "\n";
    $descendant->{Sequence} = $model->{residues}->[Samplers::weighted_random_choice($probs)];
    
    #~ print "Chosen amino acid: $descendant->{Sequence}\n";
    if (not $descendant->{HiddenState}->member($descendant->{Sequence})){
        warn "Residue is not high fitness: \n".
        "Ancestor: tree $ancestor->{Name}: $ancestor->{Sequence} " . $ancestor->{HiddenState}->as_string(). "\n".
        "Descendant: $descendant->{Name}: $descendant->{Sequence} " . $descendant->{HiddenState}->as_string()."\n";
    }
}
