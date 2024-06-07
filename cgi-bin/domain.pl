#!/usr/bin/perl

#This programme compare two files one containing id#-/+1# aac svm score and another containing domain accession no. id.
#The output file are like this id#-/+1#,aac,svm score and domain in one line.
#USE:perl domain_new.pl file1(svm) file2(domain)

open(FH,"$ARGV[0]") or die "$!";
while($line=<FH>)
{
    chomp($line);
    print "$line#";
    @array2=split('#',$line);
    #print "$array2[0] ---\n";
    open(FH1,"$ARGV[1]") or die "$!";
    while($line1=<FH1>)
    {
	chomp($line1);
	if($line1=~ m/$array2[0]/)
	    {
		@array3=split/ /,$line1;
		print "$array3[1]#";
	    }
    }
    print "\n";
    close FH1;
}
close FH;
