#!/usr/bin/perl
open(FH,"$ARGV[0]") or die "$!";
while($line=<FH>)
{
    chomp($line);
    #print "$line\n";
    @array=split(/ /,$line);
    #print "$array[3]\n";
    if(1e-5 >= $array[3])
    {
	#print "if(1e-5 >= $array[3])\n";
	print "$line\n";
    }
}
