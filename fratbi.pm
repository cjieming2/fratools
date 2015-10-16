package fratbi;

use strict;
use fralib;
use File::Basename;
use Pod::Usage;
use POSIX qw{floor ceil};

=head1 NAME

fannotatestrands

=head1 SYNOPSIS

 fannotatestrands [options] -g <2bit-file> <mk-file>

  -h        help
  -s        Alignment Score Cutoff (default: 0.8)
  -d        Alignment Score Delta Cutoff (default: 0.3)
  -g        2bit encoding of a genome
  -p	    print summary of 2bit file
  -P        print detailed summary of 2bit file
  mk-file   mk file
            a)snp-id
            b)chromosome
            c)position
            d)flanks
 
 example: fannotatestrand pscalare.mk -g pscalare.2bit
          
 Annotates the strands of flanking sequence based on chromosomal position.
 Output file is strand-annotated-<mk-file>
 
=head1 DESCRIPTION

=cut

sub connect
{	
    shift();
	my $self  = {};
	
	my $twoBitFile = shift();
	my %CHROM;
	my %META;
    my $twoBitFH;
    
	open($twoBitFH, $twoBitFile) || die "Cannot open $twoBitFile";
	
	my $value;
	my $bytesRead;
	
	$bytesRead = read($twoBitFH, $value, 4);
	my $signature = unpack('L', $value);
	$bytesRead = read($twoBitFH, $value, 4);
	my $version = unpack('L', $value);
	$bytesRead = read($twoBitFH, $value, 4);
	my $sequenceCount = unpack('L', $value);
	$bytesRead = read($twoBitFH, $value, 4);
	my $reserved = unpack('L', $value);
	
	$META{TWOBITFILE} = $twoBitFile;
	$META{FH} = $twoBitFH;
	$META{SIGNATURE} = $signature;
	$META{VERSION} = $version;
	$META{SEQUENCE_COUNT} = $sequenceCount;
	$META{RESERVED} = $reserved;
	
	for (1 .. $sequenceCount)
	{
		$bytesRead = read($twoBitFH, $value, 1);
		my $nameLength = vec($value,0,8);
		$bytesRead = read($twoBitFH, $value, $nameLength);
		my $name = unpack("A$nameLength",$value);
		$bytesRead = read($twoBitFH, $value, 4);
		my $offset = unpack('L', $value);
		
		#>gi|17981852|ref|NC_001807.4| Homo sapiens mitochondrion, complete genome
	
		if($name=~/NC_0+(\d+)/)
		{
			if($1==23)
			{
				$CHROM{X}{OFFSET} = $offset;
				$CHROM{X}{NAME} = $name;
			}
			elsif ($1==24)
			{
				$CHROM{Y}{OFFSET} = $offset;
				$CHROM{Y}{NAME} = $name;
			}
			elsif ($1==1807)
			{
				$CHROM{M}{OFFSET} = $offset;
				$CHROM{M}{NAME} = $name;
			}
			elsif ($1<=22 && $1>=1)
			{
				$CHROM{$1}{OFFSET} = $offset;
				$CHROM{$1}{NAME} = $name;
			}
			else
			{
				die "Unidentified chromosome: $1";
			}
		}
	}
	
	for my $chrom (sort(keys(%CHROM)))
	{
		seek($twoBitFH, $CHROM{$chrom}{OFFSET}, 0);
		
		$bytesRead = read($twoBitFH, $value, 4);
		my $dnaSize = unpack('L', $value);
		$bytesRead = read($twoBitFH, $value, 4);
		my $unknownBlockCount = unpack('L', $value);
		
		my @unknownBlockStart = ();
		my @unknownBlockSize = ();
		my @maskedBlockStart = ();
		my @maskedBlockSize = ();
				
		for (1 .. $unknownBlockCount)
		{
			$bytesRead = read($twoBitFH, $value, 4);
			push(@unknownBlockStart, unpack('L', $value));
		}
		
		for (1 .. $unknownBlockCount)
		{
			$bytesRead = read($twoBitFH, $value, 4);
			push(@unknownBlockSize, unpack('L', $value));
		}
		
		$bytesRead = read($twoBitFH, $value, 4);
		my $maskedBlockCount = unpack('L', $value);
		
		for (1 .. $maskedBlockCount)
		{
			$bytesRead = read($twoBitFH, $value, 4);
			push(@maskedBlockStart, unpack('L', $value));
		}
		
		for (1 .. $maskedBlockCount)
		{
			$bytesRead = read($twoBitFH, $value, 4);
			push(@maskedBlockSize, unpack('L', $value));
		}
	
		$bytesRead = read($twoBitFH, $value, 4);
		my $reservedWord = unpack('L', $value);
	
		$CHROM{$chrom}{SIZE} = $dnaSize;
		$CHROM{$chrom}{UNKNOWN_BLOCK_NO} = $unknownBlockCount;
		$CHROM{$chrom}{UNKNOWN_START} = \@unknownBlockStart;
		$CHROM{$chrom}{UNKNOWN_SIZE} = \@unknownBlockSize;
		$CHROM{$chrom}{MASKED_BLOCK_NO} = $maskedBlockCount;
		$CHROM{$chrom}{MASKED_START} = \@maskedBlockStart;
		$CHROM{$chrom}{MASKED_SIZE} = \@maskedBlockSize;
		$CHROM{$chrom}{RESERVED_WORD} = $reservedWord;			
	
		#combine blocks in order
		my $combinedBlockCount = $unknownBlockCount + $maskedBlockCount;
		my @combinedBlockStart;
		my @combinedBlockSize;
		map {push(@combinedBlockStart, $_)} @unknownBlockStart;
		map {push(@combinedBlockSize, $_)} @unknownBlockSize;
		
		for my $i (0 .. $#maskedBlockStart)
		{
			for my $j (0 .. $#combinedBlockStart)
			{
				if ($maskedBlockStart[$i]<$combinedBlockStart[$j])
				{
					splice(@combinedBlockStart, $j, 0, $maskedBlockStart[$i]);
					splice(@combinedBlockSize, $j, 0, $maskedBlockSize[$i]);
					last;
				}
				elsif ($j==$#combinedBlockStart)
				{
					push(@combinedBlockStart, $maskedBlockStart[$i]);
					push(@combinedBlockSize, $maskedBlockSize[$i]);
				}
			}
		}
		
		$CHROM{$chrom}{COMBINED_BLOCK_NO} = $combinedBlockCount;
		$CHROM{$chrom}{COMBINED_START} = \@combinedBlockStart;
		$CHROM{$chrom}{COMBINED_SIZE} = \@combinedBlockSize;	
	}

	$self->{CHROM} = \%CHROM;
	$self->{META} = \%META;
	
	bless ($self);

	return $self;
}

sub printSummary
{
	my $self = shift();
	
	my %META = %{$self->{META}};
    my %CHROM = %{$self->{CHROM}}; 
	
	print "Summary for $META{FH}\n\n";
	printf "Signature: %x\n", $META{SIGNATURE};
	printf "Version: %d\n", $META{VERSION};
	printf "Sequence Count: %d\n", $META{SEQUENCE_COUNT};
	printf "Reserved: %d\n\n", $META{RESERVED};
	
	for my $chrom (sort {if ("$a$b"=~/\D/) {$a cmp $b} else {$a <=> $b}} keys(%CHROM))
	{
	    print <<REPORT;
	=============
	Chromosome $chrom
	=============
	Name:   $CHROM{$chrom}{NAME}
	Offset: $CHROM{$chrom}{OFFSET}
	Size:   $CHROM{$chrom}{SIZE}
		
REPORT
	}		
}

sub printDetailedSummary
{
	my $self = shift();
	
	#arguments
	my $combined = shift();
	
	my %META = %{$self->{META}};
    my %CHROM = %{$self->{CHROM}}; 
	
	#prints out results
	print "Summary for  $META{FH}\n\n";
	printf "Signature: %x\n", $META{SIGNATURE};
	printf "Version: %d\n", $META{VERSION};
	printf "Sequence Count: %d\n", $META{SEQUENCE_COUNT};
	printf "Reserved: %d\n\n", $META{RESERVED};
	
	for my $chrom (sort {if ("$a$b"=~/\D/) {$a cmp $b} else {$a <=> $b}} keys(%CHROM))
	{
	    print <<REPORT;
	=============
	Chromosome $chrom
	=============
	Name:   $CHROM{$chrom}{NAME}
	Offset: $CHROM{$chrom}{OFFSET}
	Size:   $CHROM{$chrom}{SIZE}
		
REPORT

		if (!$combined)
		{
			print "Unknown Blocks ($CHROM{$chrom}{UNKNOWN_BLOCK_NO})\n";
			for my $i (1 .. $CHROM{$chrom}{UNKNOWN_BLOCK_NO})
			{
				printf "%2d) %9d %9d\n", $i, ${$CHROM{$chrom}{UNKNOWN_START}}[$i-1], ${$CHROM{$chrom}{UNKNOWN_SIZE}}[$i-1];
			}
			
			print "\nMasked Blocks ($CHROM{$chrom}{MASKED_BLOCK_NO})\n";
			for my $i (1 .. $CHROM{$chrom}{MASKED_BLOCK_NO})
			{
				printf "%2d) %9d %9d\n", $i, ${$CHROM{$chrom}{MASKED_START}}[$i-1], ${$CHROM{$chrom}{MASKED_SIZE}}[$i-1];
			}
		}
		else
		{
			print "Unknown Blocks ($CHROM{$chrom}{UNKNOWN_BLOCK_NO})\n";
			print "Masked Blocks ($CHROM{$chrom}{MASKED_BLOCK_NO})\n";
			print "Combined Blocks ($CHROM{$chrom}{COMBINED_BLOCK_NO})\n";
			for my $i (1 .. $CHROM{$chrom}{COMBINED_BLOCK_NO})
			{
				printf "%2d) %9d %9d\n", $i, ${$CHROM{$chrom}{COMBINED_START}}[$i-1], ${$CHROM{$chrom}{COMBINED_SIZE}}[$i-1];
			}
		}		
	}
}

sub getSequence
{
    my $self = shift;
    my $chromosome = shift;
    my $position = shift;
    my %CHROM = %{$self->{CHROM}};
    my %META = %{$self->{META}};
    my $twoBitFH = $META{FH};
    my $twoBitFile= $META{TWOBITFILE};

    #calculate length of sequence to extract for comparison
	my $fivePrimeLength = 60;
	my $threePrimeLength = 60;
	
	#zero based coordinates
	my $queryStart = max($position - $fivePrimeLength, 1) - 1;
	my $queryEnd = $queryStart + $fivePrimeLength + $threePrimeLength;

	my $offset = $CHROM{$chromosome}{OFFSET} + 16 + 
	            ($CHROM{$chromosome}{UNKNOWN_BLOCK_NO}+$CHROM{$chromosome}{MASKED_BLOCK_NO}) * 8;
	
	my $readStart = floor($queryStart/4);
	my $extraFivePrimeBases = $queryStart%4;
	#my $readLength = ceil(($queryEnd-$queryStart+1)/4);
	my $readLength = ((ceil($queryEnd/4)*4) - (floor($queryStart/4)*4)) / 4;
	my $extraThreePrimeBases = 4*$readLength - ($queryEnd-$queryStart+1) - $extraFivePrimeBases;
	my $snpBasePosition = $fivePrimeLength + $extraFivePrimeBases + 1;
	
	my $sequence;
	my $extractedFlanks = '';
	seek($twoBitFH, $offset + $readStart, 0) || die "Cannot seek in $twoBitFile";
	read($twoBitFH, $sequence, $readLength) || die "Cannot read $twoBitFile";
	
	my $baseNo = 0;
	
	for my $b (0 .. $readLength-1)
	{
		#silly little endian convention by !@#$ vec()
		for my $i (0 .. 3)
		{
			++$baseNo;
			if($baseNo>$extraFivePrimeBases && $baseNo<=$readLength*4-$extraThreePrimeBases)
			{
				my $base = translate(vec($sequence, $b*4+(3-$i), 2));
				
				if($baseNo == $snpBasePosition)
				{
					$extractedFlanks .= "[$base/ ]";
				}
				else
				{
					$extractedFlanks .= $base;
				}
			}
		}
	}
		
	#fill up known and masked sequence
	for my $i (0 .. $CHROM{$chromosome}{COMBINED_BLOCK_NO}-1)
	{
		my $blockStart = ${$CHROM{$chromosome}{COMBINED_START}}[$i]-1;
		my $blockEnd = $blockStart + ${$CHROM{$chromosome}{COMBINED_SIZE}}[$i] - 1;
		
		if ($blockEnd>=$queryStart && $blockStart<=$queryEnd)
		{
		    print STDERR "$extractedFlanks\n";
			warn "time to write code to fill up unknown/masked portions of the query sequence - thanks... : $chromosome:$queryStart-$queryEnd";
		}
	}
	
	return $extractedFlanks;
}

sub translate
{
	my $base = shift;
	
	if ($base==0x0)
	{
		return 'T';
	}
	elsif ($base==0x1)
	{
		return 'C';
	}
	elsif ($base==0x2)
	{
		return 'A';
	}
	elsif ($base==0x3)
	{
		return 'G';
	}	
}

#checks similarity of flanks
# <========[variation]========>
sub getFlanksSimilarity
{
	my ($flank1, $flank2) = @_;
	
	my @oneFivePrime;
	my @oneThreePrime;

	my @twoFivePrime;
	my @twoThreePrime;
		
	if ($flank1=~/([ACGTMRWSYKVHDBNX]*)\[.+\]([ACGTMRWSYKVHDBNX]*)/i)
	{
		@oneFivePrime = reverse(split(//, $1));
		@oneThreePrime = split(//, $2);
	}

	if ($flank2=~/([ACGTMRWSYKVHDBNX]*)\[.+\]([ACGTMRWSYKVHDBNX]*)/i)
	{
		@twoFivePrime = reverse(split(//, $1));
		@twoThreePrime = split(//, $2);
	}

	my $totalFivePrimeBases = min($#oneFivePrime, $#twoFivePrime)+ 1;
	my $fivePrimeConcordance = $totalFivePrimeBases;

	#compare five prime
	for my $i (0 .. $totalFivePrimeBases - 1)
	{
		if ($oneFivePrime[$i] ne $twoFivePrime[$i])
		{
			if (!baseMatch($oneFivePrime[$i],$twoFivePrime[$i]))
			{
				--$fivePrimeConcordance;
			}
		}
	}
	
	my $totalThreePrimeBases = min($#oneThreePrime, $#twoThreePrime) + 1;
	my $threePrimeConcordance = $totalThreePrimeBases;
					
	#compare three prime
	for my $i (0 .. $totalThreePrimeBases - 1)
	{
		if ($oneThreePrime[$i] ne $twoThreePrime[$i])
		{
			if (!baseMatch($oneThreePrime[$i],$twoThreePrime[$i]))
			{
				--$threePrimeConcordance;
			}
		}
	}
	
	return $totalFivePrimeBases + $totalThreePrimeBases == 0 ? 0 : ($fivePrimeConcordance + $threePrimeConcordance) / ($totalFivePrimeBases + $totalThreePrimeBases);
}

#checks similarity of flanks via alignments
# <========[variation]========>
sub getAlignedFlanksSimilarity
{
	my ($flank1, $flank2) = @_;
	
	my @oneFivePrime;
	my @oneThreePrime;

	my @twoFivePrime;
	my @twoThreePrime;

    my $oneSNP;
    my $twoSNP;

	if ($flank1=~/([ACGTMRWSYKVHDBNX]*)(\[.+\])([ACGTMRWSYKVHDBNX]*)/i)
	{
		@oneFivePrime = reverse(split(//, $1));
		$oneSNP = $2;
		@oneThreePrime = split(//, $3);
	}

	if ($flank2=~/([ACGTMRWSYKVHDBNX]*)(\[.+\])([ACGTMRWSYKVHDBNX]*)/i)
	{
		@twoFivePrime = reverse(split(//, $1));
		$twoSNP = $2;
		@twoThreePrime = split(//, $3);
	}
	
	my $totalFivePrimeBases = min($#oneFivePrime, $#twoFivePrime);
    my ($oneFivePrimeAlignment, $twoFivePrimeAlignment) = getGlobalAlignment(substr(join("", @oneFivePrime),0,scalar(@oneFivePrime)), substr(join("", @twoFivePrime),0,scalar(@twoFivePrime)));

	my $totalThreePrimeBases = min($#oneThreePrime, $#twoThreePrime);
    my ($oneThreePrimeAlignment, $twoThreePrimeAlignment) = getGlobalAlignment(substr(join("", @oneThreePrime),0,scalar(@oneThreePrime)), substr(join("", @twoThreePrime),0,scalar(@twoThreePrime)));

	my @oneFivePrimeAlignment = split(//, $oneFivePrimeAlignment);
	my @twoFivePrimeAlignment = split(//, $twoFivePrimeAlignment);
	my @oneThreePrimeAlignment = split(//, $oneThreePrimeAlignment);
	my @twoThreePrimeAlignment = split(//, $twoThreePrimeAlignment);
	
    my $total = 0;
    my $concordance = 0;
    my $gapCount = 0;
    for my $i (0 .. $#oneFivePrimeAlignment)
    {
        if ($oneFivePrimeAlignment[$i] ne '-' && $twoFivePrimeAlignment[$i] ne '-')
        {
            ++$concordance if ($oneFivePrimeAlignment[$i] eq $twoFivePrimeAlignment[$i]);
            ++$total;
        }
        else
        {
            ++$gapCount;
        }        
    }

    for my $i (0 .. $#oneThreePrimeAlignment)
    {
        if ($oneThreePrimeAlignment[$i] ne '-' && $twoThreePrimeAlignment[$i] ne '-')
        {
            ++$concordance if ($oneThreePrimeAlignment[$i] eq $twoThreePrimeAlignment[$i]);
            ++$total;
        }
        else
        {
            ++$gapCount;
        }
    }
    
    my $similarity = $total == 0 ? 0 : $concordance/$total;
    
    my $alignedSequence1 = join("", reverse(@oneFivePrimeAlignment)) . $oneSNP . join("", @oneThreePrimeAlignment);
    my $alignedSequence2 = join("", reverse(@twoFivePrimeAlignment)) . $twoSNP . join("", @twoThreePrimeAlignment);
        
    return ($similarity, $alignedSequence1, $alignedSequence2);
}

#align 2 sequences
sub getGlobalAlignment
{
    my ($sequence1, $sequence2) = @_;
    
    my @seq1 = split(//, $sequence1);
    my @seq2 = split(//, $sequence2);
    
    my @scoreMatrix = ();
    my @traceMatrix = ();
    
    #initialize borders
    $scoreMatrix[0][0] = 0;  
    $traceMatrix[0][0] = 'E';
    for my $i (1 .. scalar(@seq1))
    {
    	$scoreMatrix[0][$i] = -$i;
    	$traceMatrix[0][$i] = 'L';
    }
    for my $j (1 .. scalar(@seq2))
    {
    	$scoreMatrix[$j][0] = -$j;
    	$traceMatrix[$j][0] = 'U';
    }
    
    #score
    for my $j (1 .. scalar(@seq2))
    {
    	for my $i (1 .. scalar(@seq1))
	    {
	    	my $delta = $seq1[$i-1] eq $seq2[$j-1] ? 1 : 0;
	    	
	    	my $diag = $scoreMatrix[$j-1][$i-1] + $delta;
	    	my $up   = $scoreMatrix[$j-1][$i] - 1;
	    	my $left = $scoreMatrix[$j][$i-1] - 1;
	    	
	    	$scoreMatrix[$j][$i] = $up;
	    	$scoreMatrix[$j][$i] = $up<$left ? $left : $up;
	    	$scoreMatrix[$j][$i] = $scoreMatrix[$j][$i]<=$diag ? $diag : $scoreMatrix[$j][$i];

	    	$traceMatrix[$j][$i] = 'U';
	    	$traceMatrix[$j][$i] = $up<$left ? 'L' : 'U';
	    	$traceMatrix[$j][$i] = $scoreMatrix[$j][$i]<=$diag ? 'D' : $traceMatrix[$j][$i];

	    	#$scoreMatrix[$j][$i] = $diag;
	    	#$scoreMatrix[$j][$i] = $diag<$up ? $up : $diag;
	    	#$scoreMatrix[$j][$i] = $scoreMatrix[$j][$i]<=$left ? $left : $scoreMatrix[$j][$i];

	    	#$traceMatrix[$j][$i] = 'D';
	    	#$traceMatrix[$j][$i] = $diag<$up ? 'U' : 'D';
	    	#$traceMatrix[$j][$i] = $scoreMatrix[$j][$i]<=$left ? 'L' : $traceMatrix[$j][$i];
	    }
    }    
    
    my ($alignment1, $alignment2) = ("", "");    
    my $i = scalar(@seq1);
    my $j = scalar(@seq2);
    
    while ($i!=0 || $j!=0)
    {
    	if ($traceMatrix[$j][$i] eq 'D')
    	{
    		$alignment1 = $seq1[$i-1] . $alignment1;
    		$alignment2 = $seq2[$j-1] . $alignment2;
    		--$i;
    		--$j;
    	}
    	elsif ($traceMatrix[$j][$i] eq 'U')
    	{
    		$alignment1 = '-' . $alignment1;
    		$alignment2 = $seq2[$j-1] . $alignment2;
    		--$j;
    	}
    	elsif ($traceMatrix[$j][$i] eq 'L')
    	{
    		$alignment1 = $seq1[$i-1] . $alignment1;
    		$alignment2 = '-' . $alignment2;
    		--$i;
    	}
    }
    
    #trace matrices
    if(0)
    {
        for my $j (0 .. scalar(@seq2))
        {
        	for my $i (0 .. scalar(@seq1))
    	    {
    	        printf "\t%3d|%1s", $scoreMatrix[$j][$i], $traceMatrix[$j][$i];
    	    }
    	    
    	    print "\n";
    	}
    }
    
       
    return ($alignment1, $alignment2);
}

#returns the reverse complement of a flank
#i.e. ACGATCAGCTAAGCTCAG[A/G]ACGTGVTGATGCGT
sub reverseComplementFlanks
{
	my $flanks = shift;
	
	my @newFlanks = ();
	
	for my $base (split(//, $flanks))
	{
		unshift(@newFlanks, flankComplement($base));
	}
	
	return join("", @newFlanks);
}

#an auxiliary function that adds complements to extra symbols found in a flanking sequence
sub flankComplement
{
	my $base = shift;
	
    if ($base eq '[')
	{
		return ']';
	}
	elsif ($base eq ']')
	{
		return '[';
	}	
	elsif ($base eq '/')
	{
		return $base;
	}
	elsif ($base eq '-')
	{
		return $base;
	}
	else
	{
		return complementBase($base);
	}
}

sub max
{
	my ($val1, $val2) = @_;
	
	return $val1>$val2 ? $val1 : $val2;
}

sub min
{
	my ($val1, $val2) = @_;
	
	return $val1>$val2 ? $val2 : $val1;
}

return 1;