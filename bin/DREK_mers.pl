#!/usr/bin/perl
#usage: script input.fas output.table windowsize(default:4) threshold(default:0.5)
#not working, user defined residues: residues(default="DREKdrek")
#scan protein fasta for high DREK ratio regions and output to table

open(IN, '<', @ARGV[0]); #open fasta input
open(OUT, '>', @ARGV[1]); #write tabular output
if (@ARGV[3] != undef) { $kmer = @ARGV[3] } else { $kmer=4 } #windowsize default=4
if (@ARGV[4] != undef) { $threshold = @ARGV[4] } else { $threshold = 0.5 } #threshold default=0.5
#if (@ARGV[5] != undef) { $residues=@ARGV[5] } else { $residues="DREKdrek" } #target residues
##this not working properly, for now just define DREKdrek pattern manually in tr/ lines below

print OUT "seqID\tseq\tthreshold\tquery_residues\tcoordinates\tmotifs\tratios\n";

$fasta = do{local $/; <IN>};
@fasta = split(/\>/, $fasta);
foreach $fasrec (@fasta) {
	unless ($fasrec eq undef) {
		%pass = undef;
		@data = split(/\n/, $fasrec, 2);
		@temp = split(/\s+/, @data[0]);
		@data[0] = @temp[0];
		@data[1] =~ s/\n|\t|\s|\*//sgi;
		#initial scan - individual windows
		@pos = (1 .. length(@data[1]) );
		for $i (@pos) {
			if (($i + $kmer - 1) > length(@data[1])) { $n = length(@data[1]) } else {$n = $i + $kmer - 1}
			$temp = substr(@data[1], ($i -1), ($n - $i + 1));
			unless ($temp eq undef) {
				$ratio = ($temp =~ tr/DREKdrek//) / ($n - $i + 1);
				if ($ratio >= $threshold) {
					$pass{$i} = $n;
					##print "$i\t$n\t$temp\t$ratio\t@data[0]\t@data[1]\n";
				}
			}
		}
		#merge overlapping windows
		@last = undef;
		if (%pass ne undef) {
			$string_coords = undef;
			$string_motifs = undef;
			$string_ratios = undef;
			foreach $i (sort {$a <=> $b} keys %pass) {
				if (@last[1] eq undef) {
					@last[0] = $i;
					@last[1] = $pass{$i};
				}
				else {
					if ($i > @last[1]) {
						$x = 0;
						$y = 0;
						while (($x + $y) < 2) {
							$temp = substr(@data[1], (@last[0] -1), (@last[1] - @last[0] + 1));
							$x = (substr($temp, 0, 1) =~ tr/DREKdrek//);
							$y = (substr($temp, -1, 1) =~ tr/DREKdrek//);
							##print "$temp\t$x\t$y\n";
							if ($x == 0) {
								@last[0]++;
							}
							if ($y == 0) {
								@last[1]--;
							}
						}
						unless (length($temp) <= 1) {
							$ratio = ($temp =~ tr/DREKdrek//) / (length($temp));
							$string_coords .= @last[0] . "," . @last[1] . ";";
							$string_motifs .= $temp . ";";
							$string_ratios .= $ratio . ";";
						}
						@last[0] = $i;
						@last[1] = $pass{$i};				
					}
					else {
						@last[1] = $pass{$i};
					}
				}
			}
			$x = 0;
			$y = 0;
			while (($x + $y) < 2) {
				$temp = substr(@data[1], (@last[0] -1), (@last[1] - @last[0] + 1));
				$x = (substr($temp, 0, 1) =~ tr/DREKdrek//);
				$y = (substr($temp, -1, 1) =~ tr/DREKdrek//);
				##print "$temp\t$x\t$y\n";
				if ($x == 0) {
					@last[0]++;
				}
				if ($y == 0) {
					@last[1]--;
				}
			}
			unless (length($temp) <= 1) {
				$ratio = ($temp =~ tr/DREKdrek//) / (length($temp));
				$string_coords .= @last[0] . "," . @last[1];
				$string_motifs .= $temp;
				$string_ratios .= $ratio;
			}
		}
		$string_coords =~ s/;$//;
        $string_motifs =~ s/;$//;
        $string_ratios =~ s/;$//;
		
		print OUT @data[0] . "\t" . @data[1] . "\t" . $threshold . "\t" . $residues . "\t" . $string_coords . "\t" . $string_motifs . "\t" . $string_ratios . "\n";
	}	
}

close IN;
close OUT;