#!/usr/bin/perl

############################### Header Information ##############################
require 'cgi.perl';
use CGI;;
$query = new CGI;
&ReadParse;
print &PrintHeader;

################################ Reads Input Data ##############################
$atom = $query->param('atom');
$file = $query->param('file');
$svm_th = $query->param('svm_th');

#################Validation Of Input Sequence Data (file upload) ###################################
if($file ne '' && $atom eq '')
{
    $file=~m/^.*(\\|\/)(.*)/; 
    while(<$file>) 
    {
	$seqfi .= $_;
    }
}
elsif($atom ne '' && $file eq ''){

    $seqfi="$atom";
}

##############ACTUAL PROCESS BEGINS FROM HERE#######################
$infut_file = "/webservers/cgi-bin/submitopred";
$ran= int(rand 10000);
$dir = "/webservers/cgidocs/mkumar/temp/Ravindra/SubMitoPred/submito$ran";
system "mkdir $dir";
system "chmod 777 $dir";
$nam = 'input.'.'fasta';
open(FP1,">$dir/input_meta.fasta");
print FP1 "$seqfi\n";
#print "$seqfi\n";
close FP1;

system "/usr/bin/tr -d '\r' <$dir/input_meta.fasta >$dir/input_fi.fasta"; #Remove meta-character
system "/usr/bin/perl $infut_file/fasta.pl $dir/input_fi.fasta |/usr/bin/head -50 >$dir/input.fasta"; #Convert two line fasta file and select only 25 sequence
system "/bin/grep -c '>' $dir/input.fasta >$dir/total_seq";
system "/bin/grep '>' $dir/input.fasta |/usr/bin/cut -d '|' -f3 |/usr/bin/cut -d ' ' -f1 >$dir/protein_id"; #Grep protein id
system "/usr/local/bin/hmmscan --domtblout $dir/hmm_out -E 1e-5 $infut_file/Pfam/Pfam-A.hmm $dir/input.fasta >/dev/null";#Run hmmscan on protein sequence to find Pfam domain
system "/usr/bin/perl $infut_file/hmm.pl $dir/hmm_out |/bin/grep -v '#' >$dir/hmm_out_domain";
system "/usr/bin/perl $infut_file/evalue.pl $dir/hmm_out_domain |/usr/bin/cut -d ' ' -f1,2,3 |/usr/bin/tr '|' '#' >$dir/hmm_out_domain_evalue";
system "/usr/bin/perl $infut_file/id_compare.pl $dir/hmm_out_domain_evalue $dir/protein_id |/usr/bin/cut -d ' ' -f1,2 |/usr/bin/sort -u >$dir/domain_id";
$total_seq=`head -1 $dir/total_seq`;chomp($total_seq);
if($total_seq ne 0)
{
open (FH,"$dir/input.fasta") or die "$!";
while($line=<FH>)
{
    chomp($line);
    if($line=~ m/^>/)
    {
	$next=<FH>;
	chomp ($next);
	$n_ter=substr($next,0,25);#N-terminal residues
	$remain=substr($next,25,-25);#Remaining residues
	$c_ter=substr($next,-25);#C-terminal residues
	open(MYFILE,">$dir/sub") or die "$!";
    }
    #print MYFILE"$line\n$seq_ter\n$line\n$seq\n";
    print MYFILE"$line\n$n_ter\n$line\n$c_ter\n$line\n$remain\n";
    system "/usr/bin/perl $infut_file/aaseqformat.pl $dir/sub +1 >$dir/comp";
    open(FILE,"$dir/comp") or die "$!";
    $c=0;
    open(SAAC,">>$dir/saac") or die "$!";
    print SAAC "+1 ";
    while($file=<FILE>)
    {
	chomp($file);
	@array=split(/\+1/,$file);
	@array1=split(/ /,$array[1]);
	for($a=1;$a<=$#array1;$a++)
	{
	    $c++;
	    #print "$c:$array1[$a] ";
	    #open(SAAC,">>$dir/saac") or die "$!";
	    print SAAC "$c:$array1[$a] ";
	}
    }
    print SAAC "\n";
    close SAAC;
}
system "/usr/local/bin/svm_classify $dir/saac $infut_file/Models/level1_model $dir/svm_score_mito >/dev/null";
system "/usr/local/bin/svm_classify $dir/saac $infut_file/Models/level2_inner_mem_model $dir/svm_score_inner_mem >/dev/null";
system "/usr/local/bin/svm_classify $dir/saac $infut_file/Models/level2_outer_mem_model $dir/svm_score_outer_mem >/dev/null";
system "/usr/local/bin/svm_classify $dir/saac $infut_file/Models/level2_inter_mem_space_model $dir/svm_score_inter_mem_space >/dev/null";
system "/usr/local/bin/svm_classify $dir/saac $infut_file/Models/level2_matrix_model $dir/svm_score_matrix >/dev/null";
system "/usr/bin/paste $dir/protein_id $dir/saac $dir/svm_score_mito $dir/svm_score_inner_mem $dir/svm_score_outer_mem $dir/svm_score_inter_mem_space $dir/svm_score_matrix |/usr/bin/tr '\t' '#' >$dir/final";
system "/usr/bin/perl $infut_file/domain.pl $dir/final $dir/domain_id >$dir/final_pred";

open(FINAL_PRED,"$dir/final_pred") or die "$!";
@array=<FINAL_PRED>;
close FINAL_PRED;
for($q=0;$q<=$#array;$q++)
{
    @dom=split(/\#/,$array[$q]);
    $a=0;$b=0;$c=0;$d=0;$e=0;$f=0;
    open(UNIQUE_DOMAIN,"$infut_file/mito_non_mito_unique") or die "$!";
    while($domainfile=<UNIQUE_DOMAIN>)
    {
	chomp($domainfile);
	@domain=split(/\#/,$domainfile);
	for($x=7;$x<=$#dom;$x++)
	{
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "NON_LOCATION_DOMAIN"))
	    {
		$a=100;
	    }
	    if(("$dom[$x]" eq "$domain[0]")&&("$domain[1]" eq "MITOCHONDRIA"))
	    {
		$b=100;
		open(SUBMITO_DOMAIN,"$infut_file/submito_unique_domain") or die "$!";
		while($submitofile=<SUBMITO_DOMAIN>)
		{
		    chomp($submitofile);
		    @submito_domain=split(/\#/,$submitofile);
		    if(("$dom[$x]" eq "$submito_domain[0]")&&("$submito_domain[1]" eq "INNER_MEM"))
		    {
			$c=100;
		    }
		    if(("$dom[$x]" eq "$submito_domain[0]")&&("$submito_domain[1]" eq "OUTER_MEM"))
		    {
			$d=100;
		    }
		    if(("$dom[$x]" eq "$submito_domain[0]")&&("$submito_domain[1]" eq "INTER_MEM_SPACE"))
		    {
			$e=100;
		    }
		    if(("$dom[$x]" eq "$submito_domain[0]")&&("$submito_domain[1]" eq "MATRIX"))
		    {
			$f=100;
		    }
		}
		close SUBMITO_DOMAIN;
	    }
	}
    }
    close UNIQUE_DOMAIN;
    
    if($a==100)
    {
	open(RESULT,">>$dir/Pfam_result") or die "$!";
	print RESULT "$dom[0]\tNON_MITOCHONDRIA\tPfam_Domain\n";
	close RESULT;
    }
    if(($b==100)&&($c==0)&&($d==0)&&($e==0)&&($f==0))
    {
	#open(RESULT,">>$dir/Pfam_result") or die "$!";
	#print RESULT "$dom[0]\tMITOCHONDRIA\tDomain_based_prediction\n";
	#close RESULT;
	if(($dom[3] >= $dom[4])&&($dom[3] >= $dom[5])&&($dom[3] >= $dom[6]))
	{
	    open(SVM_RE,">>$dir/svm_result") or die "$!";
	    print SVM_RE "$dom[0]\tMITOCHONDRIAL_INNER_MEMBRANE\tSVM_Score\n";
	    close SVM_RE;
	}
	if(($dom[4] >= $dom[3])&&($dom[4] >= $dom[5])&&($dom[4] >= $dom[6]))
	{
	    open(SVM_RE,">>$dir/svm_result") or die "$!";
	    print SVM_RE "$dom[0]\tMITOCHONDRIAL_OUTER_MEMBRANE\tSVM_Score\n";
	    close SVM_RE;
	}
	if(($dom[5] >= $dom[3])&&($dom[5] >= $dom[4])&&($dom[5] >= $dom[6]))
	{
	    open(SVM_RE,">>$dir/svm_result") or die "$!";
	    print SVM_RE "$dom[0]\tMITOCHONDRIAL_INTER_MEMBRANE_SPACE\tSVM_Score\n";
	    close SVM_RE;
	}
	if(($dom[6] >= $dom[3])&&($dom[6] >= $dom[4])&&($dom[6] >= $dom[5]))
	{
	    open(SVM_RE,">>$dir/svm_result") or die "$!";
	    print SVM_RE "$dom[0]\tMITOCHONDRIAL_MATRIX\tSVM_Score\n";
	    close SVM_RE;
	}
    }
    if((($b==100))&&(($c==100)||($d==100)||($e==100)||($f==100)))
    {
	if($c==100)
	{
	    open(RESULT,">>$dir/Pfam_result") or die "$!";
	    print RESULT "$dom[0]\tMITOCHONDRIAL_INNER_MEMBRANE\tPfam_Domain\n";
	    close RESULT;
	}
	if($d==100)
	{
	    open(RESULT,">>$dir/Pfam_result") or die "$!";
	    print RESULT "$dom[0]\tMITOCHONDRIAL_OUTER_MEMBRANE\tPfam_Domain\n";
	    close RESULT;
	}
	if($e==100)
	{
	    open(RESULT,">>$dir/Pfam_result") or die "$!";
	    print RESULT "$dom[0]\tMITOCHONDRIAL_INTER_MEMBRANE_SPACE\tPfam_Domain\n";
	    close RESULT;
	}
	if($f==100)
	{
	    open(RESULT,">>$dir/Pfam_result") or die "$!";
	    print RESULT "$dom[0]\tMITOCHONDRIAL_MATRIX\tPfam_Domain\n";
	    close RESULT;
	}
    }
    if(($a==0)&&($b==0)&&($c==0)&&($d==0)&&($e==0)&&($f==0))
    {
	if($dom[2] < $svm_th)
	{
	    open(SVM_RE,">>$dir/svm_result") or die "$!";
	    print SVM_RE "$dom[0]\tNON_MITOCHONDRIA\tSVM_Score\n";
	    close SVM_RE;
	}
	if($dom[2] >= $svm_th)
	{
	    #open(SVM_RE,">>$dir/svm_result") or die "$!";
	    #print SVM_RE "$dom[0]\tMITOCHONDRIA\tSVM_based_predection\n";
	    #close SVM_RE;
	    if(($dom[3] >= $dom[4])&&($dom[3] >= $dom[5])&&($dom[3] >= $dom[6]))
	    {
		open(SVM_RE,">>$dir/svm_result") or die "$!";
		print SVM_RE "$dom[0]\tMITOCHONDRIAL_INNER_MEMBRANE\tSVM_Score\n";
		close SVM_RE;
	    }
	    if(($dom[4] >= $dom[3])&&($dom[4] >= $dom[5])&&($dom[4] >= $dom[6]))
	    {
		open(SVM_RE,">>$dir/svm_result") or die "$!";
		print SVM_RE "$dom[0]\tMITOCHONDRIAL_OUTER_MEMBRANE\tSVM_Score\n";
		close SVM_RE;
	    }
	    if(($dom[5] >= $dom[3])&&($dom[5] >= $dom[4])&&($dom[5] >= $dom[6]))
	    {
		open(SVM_RE,">>$dir/svm_result") or die "$!";
		print SVM_RE "$dom[0]\tMITOCHONDRIAL_INTER_MEMBRANE_SPACE\tSVM_Score\n";
		close SVM_RE;
	    }
	    if(($dom[6] >= $dom[3])&&($dom[6] >= $dom[4])&&($dom[6] >= $dom[5]))
	    {
	    open(SVM_RE,">>$dir/svm_result") or die "$!";
            print SVM_RE "$dom[0]\tMITOCHONDRIAL_MATRIX\tSVM_Score\n";
            close SVM_RE;
	    }
	}
    }
}
system "/bin/cat $dir/Pfam_result $dir/svm_result >$dir/Prediction";
print  "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n";
print  "<html><HEAD>\n";
print  "<TITLE>SubMitoPred::Prediction Result</TITLE>\n";
print  "<META NAME=\"description\" CONTENT=\"SubMitoPred, University of Delhi South Campus, INDIA\">\n";
print  "</HEAD><body bgcolor=\"\#FFFFE0\">\n";
print  "<h2 ALIGN = \"CENTER\"> SubMitoPred Prediction Result</h2>\n";
print  "<HR ALIGN =\"CENTER\"> </HR>\n";
print  "<p align=\"center\"><font size=4 color=black><b>The submitted protein/proteins belongs to <font color='red'></p>";
print "<table border='1' width='400' align='center'><tr><th>Protein ID</th><th>Prediction</th><th>On the Basis of</th></tr>";
}
if($total_seq == 0)
{
    system "/bin/cat $dir/Pfam_result $dir/svm_result >$dir/Prediction";
    print  "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n";
    print  "<html><HEAD>\n";
    print  "<TITLE>SubMitoPred::Prediction Result</TITLE>\n";
    print  "<META NAME=\"description\" CONTENT=\"SubMitoPred, University of Delhi South Campus, INDIA\">\n";
    print  "</HEAD><body bgcolor=\"\#FFFFE0\">\n";
    print  "<h2 ALIGN = \"CENTER\"> SubMitoPred Prediction Result</h2>\n";
    print  "<HR ALIGN =\"CENTER\"> </HR>\n";
    print  "<p align=\"center\"><font size=4 color=black><b>Please submit your sequence in fasta format</b></p>";
}
open(PRED,"$dir/Prediction") or die "$!";
while($pre=<PRED>)
{
    chomp($pre);
    @pred=split(/\t/,$pre);
    print "<tr align='center'><td>$pred[0]</td><td>$pred[1]</td><td>$pred[2]</td></tr>";
}
print "</table>";
print "</font></b></font></p>\n";
print  "<p align=\"center\"><font size=3 color=black><b>Thanks for using SubMitoPred Prediction Server</b></font></p>\n";
print  "<p align=\"center\"><font size=3 color=black><b>If you have any problem or suggestions please contact <a href='mailto:manish@south.du.ac.in'>Dr. Manish Kumar</a></b></font>. Please mention your job number in any communication.</p></br>\n";
print  "<p ALIGN=\"CENTER\"><b>Your job number is <font color=\"red\">$ran</b></font></p>\n";
print  "</body>\n";
print  "</html>\n";
system "chmod 000 $dir";
