#!/usr/bin/perl
#Usage: perl fasta-get-markov.pl -m 3 <bg_promoter.fasta >bg_promoter.model

# $Id: fasta-get-markov.txt 1339 2006-09-21 19:46:28Z tbailey $
# $Log$
# Revision 1.1  2005/07/28 23:53:19  nadya
# Initial revision
#
#
# AUTHOR: Timothy L. Bailey
# CREATE DATE: 7-1-2002

# Copyright Timothy L. Bailey, 2002

$PGM = $0;			# name of program
$PGM =~ s#.*/##;                # remove part up to last slash
@args = @ARGV;			# arguments to program
$| = 1;				# flush after all prints
$SIG{'INT'} = 'cleanup';	# interrupt handler
# Note: so that interrupts work, always use for system calls:
# 	if ($status = system($command)) {cleanup($status)}

# requires
push(@INC, split(":", $ENV{'PATH'}));	# look in entire path
#require 'beta.pl';

# defaults
$order = 0;			# order of Markov model to assume
$DNA_ALPHABET = "ACGT";		# DNA sequence alphabet
$PROT_ALPHABET = "ACDEFGHIKLMNPQRSTVWY";	# protein sequence alphabet
$ALPHABET = $DNA_ALPHABET;	# sequence alphabet
$USE_RC = 1;			# combine tuples and their reverse complements

$usage = <<USAGE;		# usage message
  USAGE:
	$PGM [-m <order>] [-p]

	[-m <order>]	order of Markov model to use; default $order
	[-p]		use protein alphabet; default: use DNA alphabet

	Estimate a Markov model from a FASTA file of sequences.

	Ignores (removes) ambiguous characters before computing
	the model.  Combines both strands of DNA.

	Reads standard input.
	Writes standard output.

        Author: Timothy L. Bailey
USAGE

$nargs = 0;					# number of required args
if ($#ARGV+1 < $nargs) { print_usage("$usage", 1); }

# get input arguments
while ($#ARGV >= 0) {
  $_ = shift;
  if ($_ eq "-h") {				# help
    print_usage("$usage", 0);
  } elsif ($_ eq "-m") {			# order of Markov model
    $order = shift;
  } elsif ($_ eq "-p") {			# use protein
    $ALPHABET = $PROT_ALPHABET;
    $USE_RC = 1;				# no reverse complements
  } else {
    print_usage("$usage", 1);
  }
}

printf stderr "Using %s alphabet...\n", ($ALPHABET eq $DNA_ALPHABET) ? "DNA" : "PROTEIN";
#
# read FASTA file
#
while (($seq = read_fasta_sequence()) ne "") {

  printf stderr "Seq: %d\r", ++$i;

  # remove the ambiguous characters
  $seq = remove_ambigs($seq, $ALPHABET);

  # count the tuples
  count_tuples($seq, $order);

} # read next sequence

# print the tuple frequencies
print_tuples($order, $ALPHABET);

# cleanup files
cleanup($status);
 
###############################################################################
#                      Subroutines                                            #
###############################################################################

###############################################################################
#	
#	read_fasta_sequence
#
#	Read a sequence from a FASTA file of sequences.
#	Removes whitespace from the sequence.
#
#	Returns the sequence (without the ID line) or "" if EOF.
#
###############################################################################
sub read_fasta_sequence {
  my($seq);

  $seq = "";					# no sequence read yet
  while (<STDIN>) {
    if (/^>/) {					# new sequence
      if ($seq ne "") { return($seq); }		# sequence already read?
    } else {
      s/\s*//g; 				# remove whitespace from string
      tr/a-z/A-Z/;				# convert to uppercase
      $seq .= $_;
    }
  }
  return($seq);					# return the sequence
} # sub read_fasta_sequence

###############################################################################
#	
#	remove_ambigs
#
#	Remove all the ambiguous characters from a sequence.
#	Works by removing all but those characters listed in the alph string.
#
#	Returns the updated sequence.
#
###############################################################################
sub remove_ambigs {
  my($seq, $alph) = @_;
  
  $seq =~ s/[^$alph]//g; 
  $seq;
} # sub remove_ambigs

###############################################################################
#
#	count_tuples
#
# 	Count the occurences of tuples of lengths 1 to order+1.
#
#	Sets global:
#	$cnt{$tuple}			counts of tuples
#	$total[$w]			total counts for width $w
#
###############################################################################
sub count_tuples {
  my($seq, $order) = @_;
  local($minw, $maxw, $length, $w, $i, $tuple);

  $minw = 1;
  $maxw = $order+1;
  $length = length($seq);

  #
  # count tuples
  #
  for ($w=$minw; $w<=$maxw; $w++) { 		# width of tuple

    for ($i=0; $i<$length-$w+1; $i++) {		# start of tuple
      $tuple = substr($seq, $i, $w);		# get the tuple
      $cnt{$tuple}++;
      $total[$w]++;
    } # start of tuple
  } # width of tuple

} # count_tuples

###############################################################################
#
#	print_tuples
#
# 	Print the frequency of tuples of lengths 1 to order+1.
#	Prints all possible tuples, including unobserved ones.
#	Combines the tuple with its reverse complement and uses average
#	frequency.
#
#	Uses global:
#	$cnt{$tuple}			counts of tuples
#	$total[$w]			total counts for width $w
#
###############################################################################
sub print_tuples {
  my($order, $alph) = @_;
  my($minw, $maxw, $w, @tuples, @new_tuples, @letters, $t, $a, $i, $tuple, $rc, $freq);
  my($sum);

  $minw = 1;
  $maxw = $order+1;
  @letters = split(//, $alph);		# letters in alphabet

  #
  # print tuple frequencies
  #
  $tuples[0] = "";
  for ($w=$minw; $w<=$maxw; $w++) { 		# width of tuple
    printf("# order %d\n", $w-1);
    $sum = 0; 
    # add each letter in alphabet to each old tuple
    $i = 0;
    foreach $t (@tuples) {
      foreach $a (@letters) {
        $tuple = $t . $a;
        if ($USE_RC) {
	  $rc = rc($tuple);
	  $freq = ($cnt{$tuple} + $cnt{$rc})/(2*$total[$w]);
        } else {
	  $freq = ($cnt{$tuple} / $total[$w]);
        }
        #$sum += $freq;
        printf("%s %9.3e\n", $tuple, $freq);
        $new_tuples[$i++] = $tuple;
      } # new letter
    } # old tuple

    @tuples = @new_tuples;

  } # width

} # print_tuples

###############################################################################
#
#	rc
#
# 	Get reverse complement of DNA string
#
###############################################################################
sub rc {
  my($string) = @_;
  my($w, $i, $seq, $first);

  $first = 0;
  if ($string =~ /^\*/) {			# handle Hamming-1 strings
    $first = 1;
    $seq = "*";			
  }
  $w = length($string);
  for ($i=$w-1; $i>=$first; $i--) {
    $a = substr($string, $i, 1);
    $a = ($a eq "A") ? "T" : ($a eq "C") ? "G" : ($a eq "G") ? "C" : 
      ($a eq "T") ? "A" : $a;
    $seq .= $a;
  }
  $seq;
} # rc

###############################################################################
#
#       print_usage
#
#	Print the usage message and exit.
#
###############################################################################
sub print_usage {
  local ($usage, $status) = @_;
 
  if (-c STDOUT) {			# standard output is a terminal
    open(C, "| more");
    print C $usage;
    close C;
  } else {				# standard output not a terminal
    print STDERR $usage;
  }

  exit $status;
}
 
###############################################################################
#       cleanup
#
#       cleanup stuff
#
###############################################################################
sub cleanup {
  local($status, $msg) = @_;
  if ($status && "$msg") {print STDERR "$msg: $status\n";}
  exit($status);
}

 
