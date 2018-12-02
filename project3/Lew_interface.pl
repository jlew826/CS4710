#storing user input as variables
$filename = $ARGV[0];
$chain1 = $ARGV[1];
$chain2 = $ARGV[2];
$inputdistance = $ARGV[3];

#opening .pdb file
if ($filename =~ /.*\.pdb/){
  open my $proteinfile, '<', $filename
    or die "Couldn't open \"$filename\"\n";

  #store all lines of the file into array @protein
  @protein=();
  while ($line = <$proteinfile>){
    chomp $line;
    $protein[@protein]=$line;
  }
  $proteinsize = @protein;

  #close the file;
  close $proteinfile;
}

else {
  print "Not a .pdb file\n";
  exit;
}

#PART 1
#Write a subroutine (findinterfaceAA) that takes as arguments the C_alpha atoms of the amino acids of each of the two
#chains and computes the distances Dist_Calpha of all pairs of C_alpha atoms, with one C_alpha in chain
#A and the other in chain B

#initializing arrays to be used
@chain1x = (); #chain 1 x-coordinates
@chain1y = (); #chain 1 y-coordinates
@chain1z = (); #chain 1 z-coordinates
@chain2x = (); #chain 2 x-coordinates
@chain2y = (); #chain 2 y-coordinates
@chain2z = (); #chain 2 z-coordinates
@chain1AAcode = (); #chain 1 AA CODE and AA NUMBER
@chain2AAcode = (); #chain 2 AA CODE and AA NUMBER

#Find all AA CODE and NUMBER of chain 1 C-alpha atms from the file and store into @chain1AAcode
$numchain1atoms = 0; #counter for the number of chain 1 C-alpha atoms
for (my $i = 0; $i < $proteinsize; $i++){
  if (@protein[$i]=~ /^ATOM.{9}CA.{6}$chain1/){
    @chain1AAcode[$numchain1atoms] = (substr @protein[$i], 21,1).":".(substr @protein[$i], 17, 3)."(".(substr @protein[$i], 23, 3).")";
    $numchain1atoms++;
  }
}

#Find all x, y, z coordinates of the chain 1 C-alpha atoms from the file,
#and store them into the array variables @chain1x, @chain1y, @chain1z, respectively
$numchain1coord = 0; #counter for the number of chain 1 coordinates
for (my $i = 0; $i < $proteinsize; $i++){
  if (@protein[$i]=~ /^ATOM.{9}CA.{6}$chain1/){
    @chain1x[$numchain1coord] = substr @protein[$i], 30, 8;
    @chain1y[$numchain1coord] = substr @protein[$i], 38, 8;
    @chain1z[$numchain1coord] = substr @protein[$i], 46, 8;
    $numchain1coord++;
  }
}

#Find all AA CODE and NUMBER of chain 2 C-alpha atms from the file and store into @chain1AAcode
$numchain2atoms = 0; #counter for the number of chain 2 C-alpha atoms
for (my $i = 0; $i < $proteinsize; $i++){
  if (@protein[$i]=~ /^ATOM.{9}CA.{6}$chain2/){
    @chain2AAcode[$numchain2atoms] = (substr @protein[$i], 21,1).":".(substr @protein[$i], 17, 3)."(".(substr @protein[$i], 23, 3).")";
    $numchain2atoms++;
  }
}

#Find all x, y, z coordinates of the chain 2 C-alpha atoms from the file,
#and store them into the array variables @chain2x, @chain2y, @chain2z, respectively
$numchain2coord = 0; #counter for the number of chain 2 coordinates
for (my $j = 0; $j < $proteinsize; $j++){
  if (@protein[$j]=~ /^ATOM.{9}CA.{6}$chain2/){
    @chain2x[$numchain2coord] = substr @protein[$j], 30, 8;
    @chain2y[$numchain2coord] = substr @protein[$j], 38, 8;
    @chain2z[$numchain2coord] = substr @protein[$j], 46, 8;
    $numchain2coord++;
  }
}

#subroutine to compute the distance between of all pairs of chain1 and chain2 C-alpha atoms
sub findinterfaceAA  {
    my ($x1_ref, $y1_ref, $z1_ref, $x2_ref, $y2_ref, $z2_ref, $chain1AAref, $chain2AAref) = @_;
  my @x1_ref = @{ $x1_ref }; #dereferencing and copying each array
  my @y1_ref = @{ $y1_ref };
  my @z1_ref = @{ $z1_ref };
  my @x2_ref = @{ $x2_ref };
  my @y2_ref = @{ $y2_ref };
  my @z2_ref = @{ $z2_ref };
  my @chain1AAref = @{ $chain1AAref };
  my @chain2AAref = @{ $chain2AAref };
  my $numchain1atoms = @x1_ref; #the total number of all chain 1 atoms
  my $numchain2atoms = @x2_ref; #the total number of all chain 2 atoms
  my $paircount = 0; #initialize pair counter

  for (my $i=0; $i < $numchain1atoms; $i++){
    for (my $j=0; $j < $numchain2atoms; $j++){
      $eucx = (@x1_ref[$i] - @x2_ref[$j])**2; #using Euclidean formula to calculate distances
      $eucy = (@y1_ref[$i] - @y2_ref[$j])**2;
      $eucz = (@z1_ref[$i] - @z2_ref[$j])**2;
      my $currentdist = sqrt($eucx + $eucy + $eucz);

      #stores chain 1 and chain 2 atom pairs in respective arrays if their distance is less than the input distance
      if ($currentdist < $inputdistance){
        @chain1pair[$paircount] = @chain1AAref[$i];
        @chain2pair[$paircount] = @chain2AAref[$j];
        $paircount++;
      }
    }
  }

#prints out the interfaced amino acids in the correct format
  for (my $i = 0; $i < $paircount; $i++){
    print "$chain1pair[$i] interacts with $chain2pair[$i]\n";
  }
  print "\n";
  return @chain1pair;
  return @chain2pair;
}

findinterfaceAA(\@chain1x, \@chain1y, \@chain1z, \@chain2x, \@chain2y, \@chain2z, \@chain1AAcode, \@chain2AAcode);

#PART 2
#Write a subroutine that for each chain computes the fraction of the interface
#amino acids lying on the secondary structures alpha helices and beta sheets.

#finding all AA on the alpha helices from .pdb file
@alphahelixlowerlimit = ();
@alphahelixupperlimit = ();
@alphahelixchain = ();
$numalphahelix = 0; #counter for the number of AA's on alphahelix
for (my $i = 0; $i < $proteinsize; $i++){
  if (@protein[$i]=~ /^HELIX/){
    @alphahelixchain[$numalphahelix] = substr @protein[$i], 19, 1;
    @alphahelixlowerlimit[$numalphahelix] = substr @protein[$i], 22, 3;
    @alphahelixupperlimit[$numalphahelix] = substr @protein[$i], 34, 3;
    $numalphahelix++;
  }
}

#finding all AA on the beta sheet from.pdb file
@betasheetlowerlimit = ();
@betasheetupperlimit = ();
@betasheetchain = ();
$numbetasheet = 0; #counter for the number of AA's on betasheet
for (my $i = 0; $i < $proteinsize; $i++){
  if (@protein[$i]=~ /^SHEET/){
    @betasheetchain[$numbetasheet] = substr @protein[$i], 21, 1;
    @betasheetlowerlimit[$numbetasheet] = substr @protein[$i], 23, 3;
    @betasheetupperlimit[$numbetasheet] = substr @protein[$i], 34, 3;
    $numbetasheet++;
  }
}

#subroutine to create a new array of only unique values from an old array.
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

@uniqchainpair1arr = uniq(@chain1pair); #only the unique values in chain 1 pair
@sortuniqchainpair1 = sort { substr($a, 6) <=> substr($b, 6)  } @uniqchainpair1arr; #sorting the array by primary seq position
@uniqchainpair2arr = uniq(@chain2pair); #only the unique values in chain 2 pair
@sortuniqchainpair2 = sort { substr($a, 6) <=> substr($b, 6)  } @uniqchainpair2arr; #sorting the array by primary seq position

#computing fraction of interface AA on alpha helices and beta sheets for chain 1
sub computefraction {
    my ($chainpair_ref) = @_;
  my @chainpair_ref = @{ $chainpair_ref }; #dereferencing and copying each array
  $chainpairlength = scalar @chainpair_ref;
  $interfacealpha = 0;
  $interfacebeta = 0;

  for (my $i =0; $i < $chainpairlength; $i++){
    for (my $j = 0; $j < $numalphahelix; $j++){
      if (((substr $chainpair_ref[$i], 0, 1) eq @alphahelixchain[$j]) &&
        ((substr $chainpair_ref[$i], 6, 3) ge @alphahelixlowerlimit[$j]) &&
        ((substr $chainpair_ref[$i], 6, 3) le @alphahelixupperlimit[$j])) {
        $interfacealpha++;
      }
    }
    for (my $j = 0; $j < $numbetasheet; $j++){
      if (((substr $chainpair_ref[$i], 0, 1) eq @betasheetchain[$j]) &&
        ((substr $chainpair_ref[$i], 6, 3) ge @betasheetlowerlimit[$j]) &&
        ((substr $chainpair_ref[$i], 6, 3) le @betasheetupperlimit[$j])) {
        $interfacebeta++;
      }
    }
  }
print (substr $chainpair_ref[0], 0, 1);
print "\n";
print "$interfacealpha/$chainpairlength of the interface amino acids lying on alpha helices.\n";
print "$interfacebeta/$chainpairlength of the interface amino acids lying on beta sheets.\n";
print "\n";
}
computefraction(\@sortuniqchainpair1);
computefraction(\@sortuniqchainpair2);

#PART 3
#Write a subroutine that for each interface amino acid of a chain determines the closest interface atom
#of the same chain to the right in the primary sequence. Then it determines the difference in position of
#the two amino acids in the sequence.

#subroutine to determine next closest IA to the right in the primary sequence
sub closestIA {
    my ($chainpair_ref) = @_;
  my @chainpair_ref = @{ $chainpair_ref }; #dereferencing and copying each array
  my $uniqchainlength = scalar @chainpair_ref;

  print (substr $chainpair_ref[0], 0, 1);
  print "\n";

  for (my $i = 0; $i < ($uniqchainlength-1); $i++) {
    my $b = substr @chainpair_ref[$i+1], 6, 3;
    my $a = substr @chainpair_ref[$i], 6, 3;
    my $distanceIA = $b - $a;
    print (substr @chainpair_ref[$i], 2, 3);
    print ": closest ".(substr @chainpair_ref[$i+1], 2, 3)." at distance $distanceIA.\n";
  }
  print "\n";
}

closestIA(\@sortuniqchainpair1);
closestIA(\@sortuniqchainpair2);
