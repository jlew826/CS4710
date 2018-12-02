#EXERCISE 1

#prompt user to enter two strings of DNA, store user input into 2 variables: DNA1 and DNA2;
#remove newline and whitespace from DNA1 and DNA2
#if string is not valid nucleotide, exit
print "Please type a (short) string of DNA: ";
$DNA1 = <STDIN>;
chomp $DNA1;
$DNA1 =~ s/\s//g;

if ($DNA1 =~ /[^acgt]/ig){
  print "This is not a valid DNA sequence.\n";
  exit;
}

print "Please type a second (short) string of DNA: ";
$DNA2 = <STDIN>;
chomp $DNA2;
$DNA2 =~ s/\s//g;

if ($DNA2 =~ /[^acgt]/ig){
  print "This is not a valid DNA sequence.\n";
  exit;
}

#count length of DNA1 (before concationation)
$lengthDNA1 = length $DNA1;

#print the two DNA strings as concatenated
$DNA1 .= $DNA2;
print ($DNA1,"\n");

#initialize index
$index = 0;

#print second string lined up over its copy at the end of the concatenated string
while ($index < $lengthDNA1){
  print " ";
  $index++;
}
print $DNA2, "\n";

#EXERCISE 2

#prompt user to enter protein sequence file name
print "Enter the protein sequence file name: ";
$proteinfileseq = <STDIN>;

# Remove the newline from the protein filename
chomp $proteinfileseq;

# open the file, or exit
unless ( open(PROTEINFILE, $proteinfileseq) ) {
    print "Cannot open file \"$proteinfileseq\"\n";
exit; }

# Read the protein sequence data from the file, and store it into the array variable @protein
@protein = ();
while ($line = <PROTEINFILE>){
  chomp $line;
  $protein[@protein]=$line;
}

# Close the file, and put all of the data from @protein into variable $protein
close PROTEINFILE;
$protein = join( '', @protein);

# Remove whitespace (if any)
$protein =~ s/\s//g;

#count length of protein
$lengthofprotein = length($protein);

#count hydrophobic amino acids (A, I, L, M, F, V, P, G)
$numhpaa = ($protein =~ tr/AILMFVPG//);

#find percentage of hydrophobic AA's
print ("The number of amino acids= ", $lengthofprotein,"\n");
print ("The number of hydrophobic amino acids= ", $numhpaa, "\n");
print ("The percentage of hydrophobic amino acids in this protein sequence is ", ($numhpaa/$lengthofprotein)*100, "%\n");


#EXERCISE 3

#prompt user to enter a DNA sequence
#remove newline and white space
#if sequence is not valid, exit
print "Enter your DNA sequence: ";
$DNAseq = <STDIN>;
chomp $DNAseq;
$DNAseq =~ s/\s//g;

if ($DNAseq =~ /[^acgt]/ig){
  print "This is not a valid DNA sequence.\n";
  exit;
}

#capitalize DNA sequence
$DNAseq = uc($DNAseq);

#find length of DNA sequence
$lenDNAseq = (length($DNAseq) * -1);

#initialize index and rcDNA
$index = -1;
$rcDNA = "";

#reads the sequence backwards and replaces nucleotide with it's complement until it gets to the beginning of DNAseq
until ($index < $lenDNAseq){
  if ((substr $DNAseq, $index, 1) eq 'A'){
    $rcDNA .= 'T';
  }
  elsif((substr $DNAseq, $index, 1) eq 'C'){
    $rcDNA .= 'G';
  }
  elsif((substr $DNAseq, $index, 1) eq 'G'){
    $rcDNA .= 'C';
  }
  elsif((substr $DNAseq, $index, 1) eq 'T'){
    $rcDNA .= 'A';
  }
  $revDNA .= substr $DNAseq, $index, 1;
  $index--;
}

print "The reverse complement of this sequence is: \n", $rcDNA, "\n";
