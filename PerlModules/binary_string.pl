use strict;
use warnings;


my $binary_string = "01000010";
print "binary string: $binary_string\n";

# Binary string to bits
my $bits = pack("B*", $binary_string);

# Bits to binary string
my $bin = unpack("B*", $bits);
print "as binary string: $bin\n";

my $hex = unpack "H*", $bits;
print "as hex: $hex\n";

my $bytes = unpack "A*", $bits;
print "as ascii chars: $bytes\n";
