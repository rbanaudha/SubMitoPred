#!/usr/bin/perl
# CONVERT MULTIPLE LINES OF FASTA FORMAT SEQUENCE INTO 2 LINE SEQUENCE
# FIRST LINE BEING THE HEADER WHILE SECOND BODY
# NOUMBER OF OUTPUT FILE =1

$f_name=$ARGV[0];

if($f_name eq '')
{
    #print "ERROR!\nfasta.pl input output\n";
    print "ERROR!\nfasta.pl input\n";
    exit;
}



open(FH,"$f_name") or print "at opening of file\n";
while($line=<FH>)
{
    chomp($line);
    $whole[$a]=$line;
    $a++;
}
close FH;

for($b=0;$b<=$#whole;$b++)
{
    if($whole[$b] =~ m/^>/)
    {
	#open(FH1,">>$out_name");
	#print FH1 "$whole[$b]\n";
	chomp($whole[$b]);
	print "$whole[$b]\n";
	$b++;
	while($whole[$b] =~ m/^[^>]/)
	{   
	    #print FH1 "$whole[$b]";
	    print "$whole[$b]";
	    $b++;
	}
	$b--;
    }
    #print FH1 "\n";
    #close FH1;
    print "\n";
}


    
