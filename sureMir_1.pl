#!/usr/bin/perl -w
#use strict;
use warnings;

use Scalar::Util qw(looks_like_number);
##use String::LCSS_XS qw(lcss lcss_all);

use Config;
    $Config{useithreads} or die('Recompile Perl with threads to run this program.');

use threads;
use Benchmark qw(:hireswallclock);

#the implementation in c of longest common subtring with modification for mismatchs using dynamic programming approach

use Inline C => q[
#include <stdio.h>
#undef free
SV * lcs_mm_sv_char_return_1(char * str1_, char * str2, int mismatch,int max_length_str1_,int min_max_length)
{
   //memory allocation for returned string
	SV * ret;
        SV * ret1;
        
	char  *str=(char *)malloc(400 * sizeof(char));
	char *str1;
        
	//we need a copy otherway arise memory acces error
	char *str_=strdup(str1_);
	
	int has_max_len=0;
	int last_has_max_len=0;
	
	
        //count in which string we have the longest alignment
	int count_max_length=0;
	int count_=0;
	
	int i = 0;
	int j = 0;
	int k = 0;
	
        //the used matrices in calcualtion
	int **num;
	int **mismatch_count;
	
	//helping variables
	int const MISMATCH_COUNT = mismatch;
	int prev_mismatch_nr = 0;
	int last_max_length = 0;
	int prev_max_length =0;
	int maxlen = 0;
	int glob_maxlen =0;
	
        //variables for counting the mismatches
	int count_mismatch=0;
        int pos_mismatch_1=0;
	
        //positions
        int maxi = 0;
	int maxj = 0;
	int max2i = 0;
	int max2j = 0;

	//this is fix now
	int len1 =  max_length_str1_;
	int len2 = strlen(str2);
	int minta_hossz = (int)(len2-2)/2;

	//memory allocation for the used matrices
	num = (int **)calloc(len1 + 1, sizeof(int *));
	for(k = 0; k < len1; k++)
		num[k] = (int *)calloc(len2 + 1, sizeof(int));
	
	mismatch_count = (int **)calloc(len1 + 1, sizeof(int *));
	for(k = 0; k < len1; k++)
        mismatch_count[k] = (int *)calloc(len2 + 1, sizeof(int));
	
        //initialization of values with emoty strings
	memset ( str, 0, 400 );
        ret1 = newSVpvn(str, 0);
        ret = newSVpvn(str, 0);
        
	//splitong of the first string (in our case mature miRNA sequences)
	str1 = strtok (str_,"*");
	
	while (str1 != NULL)
	{
   		int len1 = strlen(str1);
   
		prev_mismatch_nr = 0;
		last_max_length = 0;
		maxlen=0;

		//check if the length one of the string  is smaller then 1
		//perl inline macros and return
		//if(len1 < 1 || len2 < 1)
		// return "0";

		//initialization to value 1 for mismatch matrix (we need max 1 mismatch this can be changed for parameter value)
		for(i = 0; i < len1; i++)
                {
		    for(j = 0; j < len2; j++)
				mismatch_count[i][j]=1;
		}
		
                //initilaization to 0 the length matrix
		for(i = 0; i < len1; i++)
		{
		    for(j = 0; j < len2; j++)
				num[i][j]=0;
		}


		//counting mismatches
		count_mismatch=0;
		pos_mismatch_1=0;

		//calcualtion
		for(i = 0; i < len1; i++)
		{

	        for(j = 0; j < len2; j++)
		    {
			    if(str1[i] != str2[j])
				{	if (MISMATCH_COUNT > 0)
					{
						
						if((i == 0) || (j == 0))
						{
							num[i][j] = 1;//
						}
						else
						//we should check the values in diagonal backward from the position  num[i - 1][j - 1]
                                                //if we have 0 value greater or equal than the value of mismatch we need to start to count
                                                //the substring num[i][j] =1
                                                {
							last_max_length = num[i - 1][j - 1];
							prev_mismatch_nr = 0;
							for(k = last_max_length ; k > 0; k--)
							{
								if (mismatch_count[i - k][j - k] == 0) 
									prev_mismatch_nr++;
							}

							if (prev_mismatch_nr >= MISMATCH_COUNT)
							{
								num[i][j] =  1;//
	    
							}
							else
							{
								num[i][j] =  num[i - 1][j - 1]+1;
							}
						}//else end
				
						mismatch_count[i][j]=0;
					}//mismatch if
					else
						num[i][j] = 0;
				}//if str1...
				else
				{
					if((i == 0) || (j == 0))
						num[i][j] = 1;
					else
						num[i][j] = 1 + num[i - 1][j - 1];
				}//end if str1...
                                
				//this shoudl be checked allways not only for alignment
                                //we should have more maximal alignment-todo                               
				//we need all position which is maximal and bigger then min_max_elgth
		
				if (num[i][j]>= min_max_length)
				{
					if (num[i][j] <= maxlen)
					{
					//we need the same position with the the greatest alignment
                                        //we add position (in our case the order of mature miRNA
                                        //from the concatenated miRNA string) to the return string
                                        sv_catpvf(ret1, "*%d" ,maxi );
					sv_catpvf(ret1, "*%d" , maxj);

					sv_catpvf(ret1, "*%d" , count_mismatch);
                                        
					//this if we have  mismatch value greater then 1 we should do in a cycle-todo
                                        sv_catpvf(ret1, "*%d" , pos_mismatch_1);
					}				
				
					maxlen = num[i][j];
					has_max_len++;
					
					//checking mismatch
					count_mismatch=0;
					pos_mismatch_1=0;
					
					//if we have the mismatch greater then 1 we should save this position in cycle
					for(k =  num[i-1][j-1] ; k>=0; k--)
					{
						if (mismatch_count[i - k][j - k] == 0)
						{
							count_mismatch+=1;
							pos_mismatch_1=i-k+1;
						}
					}
					
					//saving position				
                                        //position of mature miRNA
                                        maxi = i-maxlen+2;
					maxj = i+1;
					//reads position
					max2i = j-maxlen+2;
					max2j = j+1;
					
				}
			
			}//for j

			
		}//for i

		//we add the order of mature miRNA and the max value if greater than 15 in our case
                if (has_max_len> last_has_max_len)
		{
		
			sv_catpvf(ret1, "*%d" ,maxi );
			sv_catpvf(ret1, "*%d" , maxj);
			
			sv_catpvf(ret1, "*%d" , count_mismatch);
			sv_catpvf(ret1, "*%d" , pos_mismatch_1);
			
                        sv_catpvf(ret1, "**%d" , count_);
		        sv_catpvf(ret1, "**%d*" , maxlen);
                        
			if (glob_maxlen<maxlen) glob_maxlen=maxlen; 
			last_has_max_len=has_max_len;
        	}
		
                count_++;
		str1 = strtok (NULL, "*");
     }//while
	
        //will be freed if the while cycle is ended
	for(i = 0; i < len1; i++)
        free(num[i]);
	free(num);

	for(i = 0; i < len1; i++)
        free(mismatch_count[i]);
        free(mismatch_count);
	
	free(str_);
	free(str1);

        //concat the two sv strings
	sv_catpvf(ret, "%d*" , glob_maxlen);
        sv_catsv(ret,ret1);
        
        //we free the str which is used for initialisation of ret string
        free(str);
	//decrease the reference to ret1 so the perl memoria managere eill be free it up.
	SvREFCNT_dec(ret1); 
	return ret;
 }
 
 void MyFree( void *p ) {
#undef free

    free( p );
}
//EOF
];

my $nr_mismatch = 1;
my $num_of_threads = 1;

my $results = 'result_sequence.csv';
my $maturefasta = 'mature.fasta';

my $readsnfafasta ='clean.fa';
my $readskfasta ='clean.fa';

#the reads and mature miRNA data files
open(MATURE_FASTA, $maturefasta) or die "Can't open mirbase file: $!";
open(READS_NFA_FASTA, $readsnfafasta) or die "Can't open mirbase file: $!";
open(READS_K_FASTA, $readskfasta) or die "Can't open mirbase file: $!";

#convert mature fasta file to csv and filter the hsa sequences
#read the solid files
my $first_line= '';
my $second_line= '';

my %line_fasta_csv = ();           

#mature miRNA will be concatenated into a string
my @all_mirna=();
my $max_length_mirna=0;

while (<MATURE_FASTA>) {
	chomp $_;
	my (@first_line) = split(' ' ,$_);
	
	my $second_line = <MATURE_FASTA>;
	my @result_line = () ;

	push @result_line,$second_line,$first_line[1],$first_line[4];
	chomp(@result_line);
	my $result = index($first_line[0], 'hsa');
	
	if (($result != -1) && (length($second_line)>=15))
	{
		$result_line[0]=~ s/U/T/g;
		
                #extending with columns starting now from 15 - 2column + 2x21
                for(my $i = 0;$i<=43;$i++)
		{
			push @result_line,0;	
		}		
		$line_fasta_csv{$result_line[0]} = \@result_line;  
	}
 }
 close(MATURE_FASTA);
 
my $length_of_mirna_line=0;

foreach my $key ( keys %line_fasta_csv )
{
	my (@result_line) = @{$line_fasta_csv{$key}};
        
        #we can put before the cycle...
        $length_of_mirna_line= scalar @result_line;
	
        #to be unique
	push @all_mirna,$result_line[0];
		
	if (length($result_line[0])>$max_length_mirna)
	{
	    $max_length_mirna = length($result_line[0]);
	}
} 
 
 my $all_mirna_=join("*",@all_mirna);
 

#converting the reads into a list and gong through search in mature miRNA now NFA
 my $read_first_line= '';
 my $read_second_line= '';

 my %reads_nfa_fasta_csv = ();            

while (<READS_NFA_FASTA>) {
	chomp $_;
	
	my (@read_first_line_) = split('-' ,$_);
	my (@read_first_line) = split(' ',$read_first_line_[0]);
        
	my $read_second_line = <READS_NFA_FASTA>;
	my @result_line = () ;

	push @result_line,$read_second_line,$read_first_line[1];
	chomp(@result_line);
	
	$reads_nfa_fasta_csv{$result_line[0]} = \@result_line;  #associate array to the name
 }
 close(READS_NFA_FASTA);

#convert the reads into a list and going through we search in maure miRna now K
 my $readk_first_line= '';
 my $readk_second_line= '';
 my %reads_k_fasta_csv = ();          
 
while (<READS_K_FASTA>) {
	chomp $_;
	my (@readk_first_line_) = split('-' ,$_);
	my (@readk_first_line) = split(' ',$readk_first_line_[0]);
       
	my $readk_second_line = <READS_K_FASTA>;
	my @result_line = () ;

	if (length($readk_second_line)>=15) {
            
	#filter out the elements shirter then 15
	push @result_line,$readk_second_line,$readk_first_line[1];
	chomp(@result_line);
	
	$reads_k_fasta_csv{$result_line[0]} = \@result_line;  #tömb hozzárendelése a névhez
	}
 }
 close(READS_K_FASTA);

#concat the mature miRNA seqeunces
my @thr =();
my $nr_files = 1;

my @result_all = ();

#comparing more files
for(my $m = 1;$m<=$nr_files;$m++)
{
#thread init
my $starttime = Benchmark->new;
my $finishtime;
my $timespent;

#split the read data array into the number of threads
my @nr_thread_list;
my $nr_elems_nfa = 0;
my @reads_nfa = ();

#this can be parametrizied or store into a list the files-todo
if ($m == 1)
{
	$nr_elems_nfa = scalar(keys %reads_nfa_fasta_csv);
	@reads_nfa = values %reads_nfa_fasta_csv;
}
else
{
	$nr_elems_nfa = scalar(keys %reads_k_fasta_csv);
	@reads_nfa = values %reads_k_fasta_csv;
}

my $elem_per_thread = int($nr_elems_nfa / $num_of_threads);
my $length_elem=$elem_per_thread ;
my $mod_elem=$nr_elems_nfa % $num_of_threads; 
my @reads_nfas;

#*******************
my @mirna_m_ = values %line_fasta_csv;
my $mirna_m_nr_ = @mirna_m_;

for(my $i = 1;$i<=$num_of_threads-1;$i++)
{
		my @reads_nfa_i = ();
		#splice get the length number of elements from the first parameter
		@reads_nfa_i = splice(@reads_nfa, 0, $length_elem);
		$reads_nfas[$i] = \@reads_nfa_i;
}

#elements of the last thread
if ($mod_elem != 0) {
	$length_elem =$nr_elems_nfa-$elem_per_thread*($num_of_threads-1);	
}
	my @reads_nfa_i = splice(@reads_nfa, 0, $length_elem);
	$reads_nfas[$num_of_threads] = \@reads_nfa_i;	


######################################xx
#we can give as parameters the number of threads, mismatch and target files (here mature miRna)
#and the fasta files with reads ex.: perl mismatch_solid -t 2 -m 1 -p mirna.fasta -f NFA.
#error handling-todo
#on the linux we should run only with 1 thread.probably c library problem.-todo
#using komodo or eclipse for easy development.

	#creating the threds giving parameters,gathering the result and print out toa  file or output screen.
	for(my $i = 1;$i<=$num_of_threads;$i++)
	{
	my @mirna_m = values %line_fasta_csv;
	my $mirna_m_nr = @mirna_m;
        
	#or create a paramlist
	my @ParamList = ($reads_nfas[$i], @mirna_m);
	
	#we will open so many file as many thread we have  to save separatly because if we save into the
        #same file should be ariseing some problem.
	open(GET_MIRNA_RESULTS_.$i, ">result_poz_".$i.".csv") or die "Can't write to results file: $!";
	
        $thr[$i] = threads->create({'context' => 'list'}, \&doOperation,@ParamList);
	}

	for(my $i = 1;$i<=$num_of_threads;$i++){
	     my @results = $thr[$i]->join();

	}

$finishtime = Benchmark->new;
$timespent = timediff($finishtime,$starttime);
print "\nDone!\nSpent ". timestr($timespent);

}#for file_nr

sub doOperation
{
	# Get the thread id. Allows each thread to be identified.
	my $id = threads->tid();
	my @reads_nfa_hits_ = ();
	my @InboundParameters = @_;
        
	#get the parameters from the list the reads and the mature miRNA data
        my @reads_nfa_= splice(@InboundParameters,0,1);
	my @reads_mirna_ = values %line_fasta_csv;
	
	my $mirna_m_nr = @reads_mirna_;
	my $reads_nr =  scalar @{$reads_nfa_[0]};
	my $nr_mismatch_l = $nr_mismatch;
	
	#go through the reads which we want to align to mature miRNA in aour case NFA and K sample data
	for(my $i = 0;$i<=$reads_nr-1;$i++)	
	{
	my (@line_read) = @{$reads_nfa_[0][$i]};
	
	my $length_read = length($line_read[0]);

	my $longest_substr = '';
	my $llongest_substr = 0;
	 
	#count how many  alignment is to mature miRna
	my $nr_hits_unique = 0;
	my @max_mismatch = ();

	#reverse complement-antisense
	my $revcom = reverse $line_read[0];
	$revcom =~ tr/ACGTacgt/TGCAtgca/;

	  #csak akkor nézi meg a kevesebb mismatchos esetet ha a maxra igaz
	  for (my $k = $nr_mismatch_l;$k<=$nr_mismatch_l;$k++)
	  {
		 
		my $length_mirna_match='null';
		my $name_mirna_match='null';
		my $lsubstr = "";
	        
                my $ret_val_ = lcs_mm_sv_char_return_1($all_mirna_,$line_read[0].'**'.$revcom,$k,$max_length_mirna,15);
                
		my $ret_val__ = substr($ret_val_,0,length($ret_val_)-1);
		
		my $line_=$line_read[0]."**".$revcom.",".$line_read[1].",".$ret_val__;
                
                $line_ =~ s/\n//g;
                $line_ =~ s/\r//g;
                
                print $line_;
                print "\n";
		
	  }#for mismatch
	}#for reads
	
	print "mirna nr $mirna_m_nr.\n";
	print "reads nr $reads_nr.\n";
	print "Thread $id done!\n";
	#close($RESULT_);
	
}#doOperation
