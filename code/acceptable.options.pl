our $options = {
    # Input
    tree_filename   =>  "../data/3.tree",
    sequences_filename => "../data/3x30.fasta",
    acceptabilities_filename => "../data/3x30.acceptabilities",
    generations    =>      10000,


    out_directory   =>  "../data/3x30/",

    time_between_prints =>  3,
    resample_acceptabilities_probability => 0.01,


    max_branch_length =>   0.05,

    # Sampling options
    sample_acceptabilities  => 0,

    sample_residue_acceptability_prior =>  0,
    prior_stdev   => 0.01,

    sample_switch_rate  => 0,
    switch_stdev  => 0.1,

    sample_sub_rate => 0,
    sub_stdev  => 0.1,


    # Model options
    # rate of substituting from unacceptable to acceptable
    q_a   =>  0.1,    
    # rate of substituting from acceptable to unacceptable
    q_u  =>   0.000001,
    # rate of substituting from unacceptable to unacceptable
    q_u_u  =>  0.000001,



    # Prior on acceptabiliy size
    priors   =>  [0.0000000001, 0.5, 0.00000005, 0.0001],
};

# Dependent options below

# Output
$options->{tree_out_filename} = $options->{out_directory} . "tree_out.newick";
$options->{sequences_out_filename} = $options->{out_directory} . "sequences_out.fasta";
$options->{acceptabilities_out_filename} = $options->{out_directory} . "out.acceptabilities";
