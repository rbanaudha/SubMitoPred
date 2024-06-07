#!/usr/bin/perl
#This programe extract domain from hmmscan output file.
open (FH,"$ARGV[0]") or die "$!";
while ($line=<FH>)
{
    chomp($line);
    #print "$line\n";
    @array=split(" ",$line);
    print"$array[0] ";
    print"$array[1] ";
    print"$array[3] ";
    print"$array[12]\n"
}
close FH;
