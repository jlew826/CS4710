#open .pdb file from command line, or die
$filename = $ARGV[0];
if ($filename =~ /.*\.pdb$/){
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

#PART 1: Centroid computation

#Find all x, y, z coordinates of the first chain C-alpha atoms from the file,
#and store them into the array variables @x, @y, @z, respectively
@x = ();
@y = ();
@z = ();
$numcaatoms = 0;   #counter for the number of first chain C-alpha atoms

for ($i = 0; $i < $proteinsize; $i++){
  if (@protein[$i]=~ /^ATOM.{9}CA/){
    @x[$numcaatoms] = substr @protein[$i], 31, 8;
    @y[$numcaatoms] = substr @protein[$i], 39, 8;
    @z[$numcaatoms] = substr @protein[$i], 47, 8;
    $numcaatoms++;
  }
  #ends forloop after the first instance of "TER", making sure that only
  #the first chain C-alpha atoms are used
  elsif (@protein[$i]=~ /^TER/){
    last;
  }
}

#subroutine to compute the centroid coordinates
sub computecentroid {
    my ($x_ref, $y_ref, $z_ref) = @_;
  my @x_ref = @{ $x_ref }; #dereferencing and copying each array
  my @y_ref = @{ $y_ref };
  my @z_ref = @{ $z_ref };

  #Find x, y, z coordinates of centroid, and store in array @centroid
  for (my $i=0; $i < $numcaatoms; $i++){
    $sumx += @x_ref[$i];
    $sumy += @y_ref[$i];
    $sumz += @z_ref[$i];
  }

  $centroidx = $sumx/$numcaatoms;
  $centroidy = $sumy/$numcaatoms;
  $centroidz = $sumz/$numcaatoms;

  @centroid = ($centroidx, $centroidy, $centroidz);
  return @centroid;
  }

computecentroid (\@x, \@y, \@z);

#PART 2: Distance computation

#subroutine to find distance between each c-Alpha atom and the centroid,
#and store in @distance array
sub computedistance {
    my ($x_ref, $y_ref, $z_ref, $centroid_ref) = @_;
  my @x_ref = @{ $x_ref }; #dereferencing and copying each array
  my @y_ref = @{ $y_ref };
  my @z_ref = @{ $z_ref };
  my @centroid_ref = @{ $centroid_ref };

#using Euclidean formula to calculate distances
  for (my $i=0; $i < $numcaatoms; $i++){
    $eucx = (@x_ref[$i] - $centroid_ref[0])**2;
    $eucy = (@y_ref[$i] - $centroid_ref[1])**2;
    $eucz = (@z_ref[$i] - $centroid_ref[2])**2;
    @distance[$i] = sqrt($eucx + $eucy + $eucz);
    }
  return @distance;
  }

computedistance (\@x, \@y, \@z, \@centroid);

#PART 3: Plot the results

#subroutine to create histogram
sub histogram {
    my ($distance_ref) = @_;
  my @distance_ref = @{ $distance_ref }; #dereferencing and copying array

  @sorteddist = sort { $a <=> $b } @distance_ref; #sorting the array from least to greatest
  $max = @sorteddist[$numcaatoms-1]; #finding the max within array
  $numbin = 5; #fixed number of bins is 5
  $width = (int($max/$numbin))+1; #determines bin width; "+1" accounts for values that are a decimal value greater than the integer
  my @hist = (0) x $numbin; #initializing bin number (5) of elements to 0

  #creating @hist array which will hold num of instances for each interval
  for (my $i=0; $i < $numcaatoms; $i++){
    if (@distance_ref[$i] < $width){
      @hist[0]++;
    }
    elsif (@distance_ref[$i] < ($width*2)){
      @hist[1]++;
    }
    elsif (@distance_ref[$i] < ($width*3)){
      @hist[2]++;
    }
    elsif (@distance_ref[$i] < ($width*4)){
      @hist[3]++;
    }
    elsif (@distance_ref[$i] < ($width*5)){
      @hist[4]++;
    }
  }

#prints the histogram along with the numerical number of distances in that range displayed after "#"'s
  for (my $j = 0; $j < $numbin; $j++){
    print (($j+1)*$width, "\t");
    for (my $k = 0; $k < @hist[$j]; $k++){
      print "#";
    }
    print (" (", "@hist[$j]", ")");
    print "\n";
  }
}
histogram (\@distance);
