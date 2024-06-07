#!/usr/bin/perl

#Extract amino acid composition from fasta file
#USAGE: perl aaseqformat.pl <file name> +/-1
#Format: 1:aaper,2:aaper....

open(FH,"$ARGV[0]") or die "$!";
while($line=<FH>)
{
    chomp($line);
    $count = 0;
    #print "$line\n";
    $array[$a]=$line;
    $a++;
}
#print "=====@array\n";
for($b=0;$b<=$#array;$b++)
{
    if($array[$b]=~ m/^>/)
    {
	#print "\n$array[$b]\n";
	#print "$ARGV[1] ";
	#$b++;
	#@array_seq1=split(/\|/,$array[$b]);
	#print ">$array_seq1[1]\t";
        #print "\n@array_seq1\n";
	@array_seq2=split(/ /,$array_seq1);
	print "$array_seq2";
    }
    while(($array[$b+1] !~ m/^>/)&&($b<=$#array))
    {
	$seq="$seq"."$array[$b+1]";
	#print "===$seq\n";
	$b++;
    }
    #print "\n$b";
    #print "++++++++$seq\n";
    @array_seq=split(//,$seq);
    #print "------@array_seq\n";
    
    for($c=0;$c<=$#array_seq;$c++)
    {
	#print "$array_seq[$c]\t";
	if($array_seq[$c] eq 'A')
	{
	    $count_A++;
	    #print "++$count_A\t";
	}
	elsif($array_seq[$c] eq 'C')
	{
	    $count_C++;
	    #print "++$count_C\t";
	}
	elsif($array_seq[$c] eq 'D')
	{
	    $count_D++;
	    #print "++$count_D\t";
	}
	elsif($array_seq[$c] eq 'E')
	{
	    $count_E++;
	    #print "++$count_E\t";
	}
	elsif($array_seq[$c] eq 'F')
	{
	    $count_F++;
	    #print "++$count_F\t";
	}
	elsif($array_seq[$c] eq 'G')
	{
	    $count_G++;
	    #print "++$count_G\t";
	}
	elsif($array_seq[$c] eq 'H')
	{
	    $count_H++;
	    #print "++$count_H\t";
	}
	elsif($array_seq[$c] eq 'I')
	{
	    $count_I++;
	    #print "++$count_I\t";
	}
	elsif($array_seq[$c] eq 'K')
	{
	    $count_K++;
	    #print "++$count_K\t";
	}
	elsif($array_seq[$c] eq 'L')
	{
	    $count_L++;
	    #print "++$count_L\t";
	}
	elsif($array_seq[$c] eq 'M')
	{
	    $count_M++;
	    #print "++$count_M\t";
	}
	elsif($array_seq[$c] eq 'N')
	{
	    $count_N++;
	    #print "++$count_N\t";
	}
	elsif($array_seq[$c] eq 'P')
	{
	    $count_P++;
	    #print "++$count_P\t";
	}
	elsif($array_seq[$c] eq 'Q')
	{
	    $count_Q++;
	    #print "++$count_Q\t";
	}
	elsif($array_seq[$c] eq 'R')
	{
	    $count_R++;
	    #print "++$count_R\t";
	}
	elsif($array_seq[$c] eq 'S')
	{
	    $count_S++;
	    #print "++$count_S\t";
	}
	elsif($array_seq[$c] eq 'T')
	{
	    $count_T++;
	    #print "++$count_T\t";
	}
	elsif($array_seq[$c] eq 'V')
	{
	    $count_V++;
	    #print "++$count_V\t";
	}
	elsif($array_seq[$c] eq 'W')
	{
	    $count_W++;
	    #print "++$count_W\t";
	}
	elsif($array_seq[$c] eq 'Y')
	{
	    $count_Y++;
	    #print "++$count_Y\t";
	}
    }
    #print "A=$count_A $#array_seq\n";
    $per_A=($count_A*100)/($#array_seq+1);
    #print "Aper=$per_A\n";
    print "$ARGV[1] ";
    printf("%4.2f ", $per_A);
    #print "2:";
    $count_A=0;
    $seq='';
    
#print "C=$count_C $#array_seq\n";
    $per_C=($count_C*100)/($#array_seq+1);
    #print "Cper=$per_C\n";
    printf("%4.2f ", $per_C);
    #print "3:";
    $count_C=0;
    $seq='';
    
#print "D=$count_D $#array_seq\n";
    $per_D=($count_D*100)/($#array_seq+1);
    #print "Dper=$per_D\n";
    printf("%4.2f ", $per_D);
    #print "4:";
    $count_D=0;
    $seq='';
    
#print "E=$count_E $#array_seq\n";
    $per_E=($count_E*100)/($#array_seq+1);
    #print "Eper=$per_E\n";
    printf("%4.2f ", $per_E);
    #print "5:";
    $count_E=0;
    $seq='';
    
#print "F=$count_F $#array_seq\n";
    $per_F=($count_F*100)/($#array_seq+1);
    #print "Fper=$per_F\n";
    printf("%4.2f ", $per_F);
    #print "6:";
    $count_F=0;
    $seq='';
#print "G=$count_G $#array_seq\n";
    $per_G=($count_G*100)/($#array_seq+1);
    #print "Gper=$per_G\n";
    printf("%4.2f ", $per_G);
    #print "7:";
    $count_G=0;
    $seq='';
    
#print "H=$count_H $#array_seq\n";
    $per_H=($count_H*100)/($#array_seq+1);
    #print "Hper=$per_H\n";
    printf("%4.2f ", $per_H);
    #print "8:";
    $count_H=0;
    $seq='';
    
#print "I=$count_I $#array_seq\n";
    $per_I=($count_I*100)/($#array_seq+1);
    #print "Iper=$per_I\n";
    printf("%4.2f ", $per_I);
    #print "9:";
    $count_I=0;
    $seq='';
    
#print "K=$count_K $#array_seq\n";
    $per_K=($count_K*100)/($#array_seq+1);
    #print "Kper=$per_K\n";
    printf("%4.2f ",$per_K);
    #print "10:"; 
    $count_K=0;
    $seq='';
    
#print "L=$count_L $#array_seq\n";
    $per_L=($count_L*100)/($#array_seq+1);
    #print "Lper=$per_L\n";
    printf("%4.2f ", $per_L);
    #print "11:";
    $count_L=0;
    $seq='';
    
#print "M=$count_M $#array_seq\n";
    $per_M=($count_M*100)/($#array_seq+1);
    #print "Mper=$per_M\n";
    printf("%4.2f ", $per_M);
    #print "12:";
    $count_M=0;
    $seq='';
    
#print "N=$count_N $#array_seq\n";
    $per_N=($count_N*100)/($#array_seq+1);
    #print "Nper=$per_N\n";
    printf("%4.2f ", $per_N);
    #print "13:";
    $count_N=0;
    $seq='';
    
#print "P=$count_P $#array_seq\n";
    $per_P=($count_P*100)/($#array_seq+1);
    #print "Pper=$per_P\n";
    printf("%4.2f ", $per_P);
    #print "14:";
    $count_P=0;
    $seq='';
    
#print "Q=$count_Q $#array_seq\n";
    $per_Q=($count_Q*100)/($#array_seq+1);
    #print "Qper=$per_Q\n";
    printf("%4.2f ", $per_Q);
    #print "15:";
    $count_Q=0;
    $seq='';
    
#print "R=$count_R $#array_seq\n";
    $per_R=($count_R*100)/($#array_seq+1);
    #print "Rper=$per_R\n";
    printf("%4.2f ", $per_R);
    #print "16:";
    $count_R=0;
    $seq='';
    
#print "S=$count_S $#array_seq\n";
    $per_S=($count_S*100)/($#array_seq+1);
    #print "Sper=$per_S\n";
    printf("%4.2f ", $per_S);
    #print "17:";
    $count_S=0;
    $seq='';
    
#print "T=$count_T $#array_seq\n";
    $per_T=($count_T*100)/($#array_seq+1);
    #print "Tper=$per_T\n";
    printf("%4.2f ", $per_T);
    #print "18:";
    $count_T=0;
    $seq='';
    
#print "V=$count_V $#array_seq\n";
    $per_V=($count_V*100)/($#array_seq+1);
    #print "Vper=$per_V\n";
    printf("%4.2f ", $per_V);
    #print "19:";
    $count_V=0;
    $seq='';
    
#print "W=$count_W $#array_seq\n";
    $per_W=($count_W*100)/($#array_seq+1);
    #print "Wper=$per_W\n";
    printf("%4.2f ", $per_W);
    #print "20:";
    $count_W=0;
    $seq='';
    
#print "Y=$count_Y $#array_seq\n";
    $per_Y=($count_Y*100)/($#array_seq+1);
    #print "Yper=$per_Y\n";
    printf("%4.2f\n", $per_Y);
    #print "$ARGV[1] ";
    #print "1:";
    $count_Y=0;
    $seq='';
}
