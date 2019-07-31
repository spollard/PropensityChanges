# Package for handling nucleotide or amino acid sequences

package Sequences;

use strict;
use warnings;


unless (caller) {
    
    my $s = FromFasta("10seqs_2.fasta");
    ToNexus($s, "10taxa.nexus");
    exit;
    
    my $sequences = FromFasta("../Data/seqs.fasta");
    print Dumper($sequences);
}

sub FromFasta {
    my ($filename) = @_;

    open my $fh, "<", $filename or die $!;

    my $sequences = {};
    
    my $name = "";
    my $sequence = "";

    while (<$fh>) {
        if (/\/\//) {
            warn "Found multiple generations of sequences in $filename"; 
            last;
        }
        chomp;
        if (s/>//) {
            if ($sequence ne "") {
                $sequences->{$name} = $sequence;
                $sequence = "";
            }
            $_ =~ s/\r//g;
            $name = $_;
        }
        else { $sequence .= $_ }
    }

    # For the last sequence
    $sequences->{$name} = $sequence;

    return $sequences;
}

sub ToFasta {
    my ($sequences, $fasta_file) = @_;
    open my $fh, ">", $fasta_file; 
    
    while (my ($sequence_name, $sequence) = each %$sequences) {
        print $fh ">", $sequence_name , "\n" , $sequence, "\n"; 
    }
}

sub ToNexus {
    my ($sequences, $nexus_file, $type) = @_;
    
    $type = $type // 'protein'; #assume protein sequence type
    
    my $number_of_taxa = scalar keys %$sequences;
    my $length_of_sequences = length $sequences->{(keys %$sequences)[0]};
    
    
    open my $fh, ">", $nexus_file; 
    print $fh "#NEXUS
    BEGIN DATA;
    dimensions ntax=$number_of_taxa nchar=$length_of_sequences;
    format datatype=$type missing=? gap=-;
    matrix\n";
    
    while (my ($sequence_name, $sequence) = each %$sequences) {
        print $fh "\t\t", $sequence_name, "\t", $sequence, "\n";    
    }
    print $fh "\n;\nEND;\n";
    
}

sub ToPhylip {
    my ($sequences, $nexus_file, $type) = @_;
    
    die "Use http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_phylip.php for this conversion";
    
    $type = $type // 'protein'; #assume protein sequence type
    
    my $number_of_taxa = scalar keys %$sequences;
    my $length_of_sequences = length $sequences->{(keys %$sequences)[0]};
    
    
    open my $fh, ">", $nexus_file; 
    print $fh "#NEXUS
    BEGIN DATA;
    dimensions ntax=$number_of_taxa nchar=$length_of_sequences;
    format datatype=$type missing=? gap=-;
    matrix\n";
    
    while (my ($sequence_name, $sequence) = each %$sequences) {
        print $fh "\t\t", $sequence_name, "\t", $sequence, "\n";    
    }
    print $fh "\n;\nEND;\n";
    
}

sub ToStockholm {
    my ($sequences, $stockholm_file) = @_;
    
    open my $fh, ">", $stockholm_file; 
    print $fh "# STOCKHOLM 1.0\n\n";
    
    while (my ($sequence_name, $sequence) = each %$sequences) {
        print $fh $sequence_name, "\t", $sequence, "\n";    
    }
    print $fh "\n//\n";
}

sub NextFromSeqoutfile {
    my ($fh) = @_;
    
    my $sequences = {};
    
    my $name = "";
    my $sequence = "";

    while (<$fh>) {
        if (/\/\//) {
            last;
        }
        chomp;
        if (s/>//) {
            if ($sequence) {
                $sequences->{$name} = $sequence;
                $sequence = "";
            }
            $name = $_;
        }
        else { $sequence .= $_ }
    }

    # For the last sequence
    $sequences->{$name} = $sequence;

    return $sequences;
}

sub FromRichards {
    my ($filename) = @_;

    open my $fh, "<", $filename or die $!;

    my $sequences = {};
    
    # kkk   628 Ichthyophis_glutinosus_194415:0.1047862 IPQRHASEKR
    while (<$fh>) {
        my $fields = [split];
        if ($fields->[0] eq "ggg" or $fields->[0] eq "kkk") {
            $sequences->{$fields->[2]} = $fields->[3];
        }
    }
    
    return $sequences;
}


sub FromRichards_LeavesOnly {
    my ($filename) = @_;

    open my $fh, "<", $filename or die $!;

    my $sequences = {};
    
    # kkk   628 Ichthyophis_glutinosus_194415:0.1047862 IPQRHASEKR
    while (<$fh>) {
        my $fields = [split];
        if ($fields->[0] eq "kkk") {
            $sequences->{$fields->[2]} = $fields->[3];
        }
    }
    
    return $sequences;
}

sub Combine {
    my ($sequences_1, $sequences_2, $sep) = @_;
    
    $sep = $sep // "";
    
    my $names_1 = [keys %$sequences_1];
    my $names_2 = [keys %$sequences_2];
    
    for my $name_1 (@$names_1) {
        if (not( $name_1 ~~ @$names_2)) {
            warn "$name_1 not found in second set of sequences";    
        }
    }
    
    for my $name_2 (@$names_2) {
        if (not( $name_2 ~~ @$names_1)) {
            warn "$name_2 not found in first set of sequences"; 
        }
        $sequences_1->{$name_2} .= $sep . $sequences_2->{$name_2}; 
    }
    
    return $sequences_1;
}

sub CombineFast {
    # Assumes that the names match. 
    my ($sequences_1, $sequences_2) = @_;
    
    for my $name (keys %$sequences_1) {
        $sequences_1->{$name} .= $sequences_2->{$name}; 
    }
    
    return $sequences_1;
}

sub Bootstrap {
    my ($sequences) = @_;
    
    my $length = length($sequences->{(keys %$sequences)[0]});
    my $sites = [map {int rand $length} 1..$length];
    
    my $bootstrap = Select($sequences, $sites);
    
    return $bootstrap, $sites
}

sub Select {
    my ($sequences, $sites) = @_;
    
    my $selection = {};
    for my $name (keys %$sequences) {
        $selection->{$name} = "";
        for my $site (@$sites) {
            if ($site <=0) {die "Selection site is <= 0"}
            if ($site > length($sequences->{$name})) {die "Selection site is beyond end of sequence"}
            $selection->{$name} .= substr $sequences->{$name}, $site-1, 1;
        }
    }
    
    return $selection
}

sub MaskFromSites {
    my ($sites) = @_;
    die "MaskFromSites is not implemented";
    my $mask;
    return $mask
    
}

sub SitesFromMask {
    my ($mask) = @_;
    
    my $sites = [];
    while ( my ($index, $value) = each @$mask) {
        if ($value) {
            push @$sites, ($index) x $value;
        }
    }
    return $sites;
}
