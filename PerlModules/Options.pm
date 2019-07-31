# Reads a control file into options and writes options to a control file.  
#
# Stephen Pollard
# 2014-8-5


package Options;

use strict;
use warnings;
use diagnostics;

use File::Copy "copy";

unless (caller) {

}

sub Read {
    my ($control_file) = @_;
    
    my $options = {};
    
    open my $fh, "<", $control_file or die $!; 
    while (<$fh>) {
        chomp;
        s/\n//g;
        s/\r//g;
        next if not $_;
        my $fields = [split " ", $_, 2]; 
        next if not scalar @$fields;
        next if $fields->[0] eq "#";
        # Some options can be arrays or hashes
        if (substr($fields->[1],0,1) ne '"' and substr($fields->[1],0,1) ne "'" 
            and ($fields->[1] =~ /\[/ or $fields->[1] =~ /\{/ )){
                $options->{$fields->[0]} = eval($fields->[1]);
            }
        else {
            $options->{$fields->[0]} = $fields->[1];
        }
    }
    
    return $options;
}

sub Write {
    my ($options, $control_file) = @_;
    
    open my $fh, ">", $control_file; 
    
    while (my ($option, $value) = each %$options) {
        print $fh "$option\t$value\n";
    }
}

sub Append {
    my ($control_file_template, $options, $control_file) = @_;
    
    copy $control_file_template, $control_file;
    
    open my $fh, ">>", $control_file; 
    print $fh "\n";
    while (my ($option, $value) = each %$options) {
        print $fh "$option\t$value\n";
    }
}

1;
