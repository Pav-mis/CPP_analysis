open(IN, '<', $ARGV[0]) or die "Can't open input file: $!";
open(OUT, '>', $ARGV[1]) or die "Can't open output file: $!";
while (<IN>) {
    $line = $_;
    chomp $line;
    unless ($line eq undef || $line =~ m/motif_coords/) {
        @data = split(/\t/, $line);
        @temp = split(/\s+/, $data[0]);
        $id = $data[0];
        $raw_seq = $data[1];
        $raw_length = length($raw_seq);
        $cutsite = 0;
        if ($data[3] > $cutsite) {
            $cutsite = $data[3];
        }
        if ($data[5] > $cutsite) {
            $cutsite = $data[5];
        }
        $mature_seq = substr($raw_seq, $cutsite);
        $mature_length = length($mature_seq);
        @cys = split /C/, $mature_seq;
        $count = ($mature_seq =~ tr/C/C/);
        $pattern = undef;
        $lengths = undef;
        $lastlen = 0;
        my $first_element = 1;
        foreach $section (@cys) {
            $len = length($section);
            if ($first_element) {
                $len += $cutsite; # Add 5 to the length of the first element
                $first_element = 0; # Update the flag
            }
            if ($lastlen > 0) {
                $pattern .= "-";
            }
            unless ($lengths eq undef) {
                $pattern .= "C";
                $lengths .= ".";
            } 
            $lengths .= $len;
            $lastlen = $len;
        }
        if ($len > 0) {
            $pattern .="-";
        }

        print OUT ">" . $id . "\n" . $mature_seq . "\n" ;
    }
}
close(IN);
close(OUT);