#!/usr/bin/perl
#This programme compare two file
#First file formate Actin PF00022.14 sp#P61162#ACTZ_CANFA_CENTROMERE
#Second file formate >sp|P61162|ACTZ_CANFA_CENTROMERE
open(FH,"$ARGV[0]") or die "$!";
while($line=<FH>)
{
    chomp($line);
    @array=split(/#/,$line);
    #print "$array[1]\n";
    open(FH1,"$ARGV[1]") or die "$!";
    while($line1=<FH1>)
    {
	chomp($line1);
	if($line1=~ m/$array[2]/)
	{
	    print "$line1 $array[0]\n";
	}
    }
}
close FH;
