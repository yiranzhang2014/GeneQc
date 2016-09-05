open(my $fh, '>', 'timereport.txt'); #time start#####################################
print $fh "Blast start\n";
$Start = time();

$geneseq = $ARGV[0];
$gff=$ARGV[1];
$bamfile = $ARGV[3];
$outfile = $ARGV[2];
mkdir "./temp/";
open MB,">./temp/log.txt";
mkdir $outfile;
use File::Basename;
use Cwd qw/abs_path/;
$Script_dir = dirname(abs_path($0));
$blast = $Script_dir."/ncbi-blast-2.2.27+-x64-linux/ncbi-blast-2.2.27+/bin/";
$samtools = $Script_dir."/samtools-1.3.1/samtools";

$makeblastdb = $blast."makeblastdb";
system("$makeblastdb -in $geneseq -out ./temp/self -dbtype nucl");
$blastn = $blast."blastn";
mkdir "./temp/blast/";
system("$blastn -query $geneseq -db ./temp/self -out ./temp/blast/self -outfmt 6 -evalue 0.000001");
print MB "$blastn -query $geneseq -db ./temp/self -out ./temp/blast/self -outfmt 6 -evalue 0.000001\n";

open CA,"./temp/blast/self";

while($line1 = <CA>)
{
        @array1 = split /\t/,$line1;
		 $array1[0] = substr ($array1[0],0,-2);
		 $array1[1] = substr ($array1[1],0,-2);
		print $array1[0]."\n"; 

        if(($array1[0] !~ $array1[1]) && ($array1[2]>90))
        {
	        if(not defined ($tier1{$array1[0]}))
                {
                        $tier1{$array1[0]}=1;
                }
                if((not defined ($group{$array1[0]}{$array1[1]})) or (not defined ($group{$array1[1]}{$array1[0]})))
                {
                        $group{$array1[0]}{$array1[1]}=1;
                }
        }
}
close CA;

open CB,">./temp/tier1.txt";
foreach $keys ( keys %tier1 )
{
        print CB $keys."\n";
}
close CB;

open CB,">./temp/group.txt";
foreach $keys (keys %group)
{
        foreach $key2 ( keys %{$group{$keys} })
        {
				
                print CB $keys."\t".$key2."\n";
        }       
}
close CB;
%tier1=();%group=();

print MB "blast done\n";

$End = time();
$Diff = $End - $Start;#############################
print $fh "Blast end. Time use:".($Diff / 3600)."h\n";

print $fh "extranct start\n";
$Start = time();

open CA,"$gff";
open CB,">./temp/gene_table.txt";
while($line1 = <CA>)
{
        @array1 = split /\s+/,$line1;
        print $array1[2]."\n";
	if($array1[2]=~/gene/)
        {
                chomp($line1);
                @array2 = split /Name=/,$array1[8];
                print CB $array2[1]."\t".$array1[0]."\t".$array1[3]."\t".$array1[4]."\n";
        }
}
close CA;
close CB;

open CA,"./temp/gene_table.txt";
while($line1 = <CA>)
{
        chomp($line1);
        @array1 = split /\s+/,$line1;
        $loc = $array1[1].":".$array1[2]."-".$array1[3];
        system("$samtools view $bamfile $loc >$outfile$array1[0]"); 
	$time =localtime;
	print MB $array1[0]."\t".$time."\n";

}
close CA;

print MB "extranct done\n";

$End = time();
$Diff = $End - $Start;#############################
print $fh "ExTranct end. Time use:".($Diff / 3600)."h\n";
print $fh "Step3 start\n";
$Start = time();

opendir(AA,$outfile);
my @list = readdir(AA);
closedir AA;
open CB,">./test/count_unique";
foreach $file(@list)
{
        if($file !~ /^\./)
        {
                open CA,"$outfile$file";
                $first_unique=$all=0;
                while($line1 = <CA>)
                {
                        @array1 = split /\s+/,$line1;
                        if($array1[1]<100)
                        {
                                $first_unique++;
                        }
                        if(not defined $gene{$array1[0]})
                        {
                                $gene{$array1[0]}=1;
                                $all++;
                        }
                }
                if($all>0)
                        {
                        $avg = $first_unique/$all;
                        }
                        else
                        {
                        $avg=1;
                        }
                print CB $file."\t".$first_unique."\t".$all."\t".$avg."\n";
                close CA;
				%gene=();
        }
}
close CB;

open CA,"./temp/tier1.txt";
while($line1 = <CA>)
{
	chomp($line1);
	@array1 =split /\s+/,$line1;
	$c1{$array1[0]}=1;
}
close CA;
open CA,"./test/count_unique";
while($line1 = <CA>)
{
	chomp($line1);
	@array1 = split /\s+/,$line1;
	$total{$array1[0]}=$array1[2];
	if($array1[3]>0.8)
	{
		$c2{$array1[0]}=1;
	}
	if($array1[3]<0.5)
	{
		$c3{$array1[0]}=1;
	}
	else {
		$c4{$array1[0]}=1;
	}
}
close CA;

open CA,"./temp/group.txt";
while($line1 = <CA>)
{
	chomp($line1);
	@array1 = split /\s+/,$line1;
	open CAA,"$outfile$array1[0]";
	while($line2 = <CCA>)
	{
		@array2 = split /\s+/,$line2;
		$count{$array1[0]}=1;
	}
	close CAA;
	open CAA,"$outfile$array1[1]";
	$number=0;
	while($line2 = <CCA>)
	{
		@array2 = split /\s+/,$line2;
		if((defined $count{$array1[0]}) && (not defined $count2{$array1[0]}))
		{
			$count2{$array1[0]}=1;
			$number++;
		}
	}
	close CAA;
	$count=();$count2=();
	if(not defined $mul{$array1[0]})
	{
		$mul{$array1[0]}=$number;
	}
	else {
		if($mul{$array1[0]}<$number)
		{
			$mul{$array1[0]}=$number;
		}
	}
	if(not defined $mul{$array1[1]})
		{
			$mul{$array1[1]}=$number;
		}
		else {
			if($mul{$array1[1]}<$number)
			{
				$mul{$array1[1]}=$number;
			}
		}
}
close CA;

foreach $var(keys %total)
{
	if(defined $c1{$var}) 
		{
			if( defined $c2{$var})
			{
				print $var."\tB\n";
			}
			if( defined $c2{$var})
			{
				$k=$total{$var}*0.5;
				if($k<$mul{$var})
				{
					print $var."\tC\n";
				}
			}
		}
	else{
		if(defined $c2{$var})
			{
				print $var."\tA\n";
			}
	
		else {
		 	print $var."\tD\n";
		}
	}
}

$End = time();
$Diff = $End - $Start;#############################
print $fh "Step3 end. Time use:".($Diff / 3600)."h\n";

close $fh;