/******************************************************************************/
/*  runny.c  */ 
/*  code to calculate homopeptide content and its significance */  
/****  Copyright 2025. Paul Martin Harrison. ****/ 
/****
 ****  Licensed under the 3-clause BSD license. See LICENSE.txt bundled with the package. 
 ****/ 
/****  
 ****  to compile: 
 ****   gcc -O2 -o runny runny.c -lm 
 ****
 ****  to run runny and get help helpy: 
 ****   ./runny -h 
 ****
****  This program is part of patterny
 ****
 ****  The latest version of this code is available at:
 ****    https://github.com/pmharrison
 **** 
 ****  Citations: 
 ****    Wang, Y, & Harrison, PM. Homopeptide and homocodon levels across fungi are 
 ****     coupled to GC/AT-bias and intrinsic disorder, with unique behaviours for 
 ****     some amino acids (2021), 11(1):10025. doi: 10.1038/s41598-021-89650-1  
 ****
 ****    Harrison, PM. "Patterny: A troupe of decipherment helpers for intrinsic disorder, 
 ****     low complexity and compositional bias in proteins", Biomolecules, ... 
 ****/ 
/*****************************************************************************************/ 

#include <stdio.h> 
#include <string.h> 
#include <stdlib.h> 
#include <math.h> 
#include <ctype.h> 
#include <sys/types.h>
#include <sys/stat.h> 
#include <time.h> 
#include <unistd.h>
#include <sys/time.h> 
#include <sys/resource.h> 
#include <stdbool.h>

#define AA_ALPHABET_SIZE 21 
#define MAXIMUM_LENGTH 5000
#define VERSION_NUMBER 1.0
#define MAX_SEQUENCE_NAME_SIZE 200 
#define MAX_FILE_NAME_SIZE 500 
#define SCRAMBLED_SAMPLE 1000 

int alphabet_length=AA_ALPHABET_SIZE, sequence_length=MAXIMUM_LENGTH, number_of_sequences; 
int *whole_rescounts, *whole_sorted_rescounts; 
int *residue_hpep, max_sequence_len, min_length=MAXIMUM_LENGTH; ; 

char *sequence, *cleanseq, *randseq, *max_sequence, sequence_name[MAX_SEQUENCE_NAME_SIZE+1] = {"initial"}; 
char line[MAXIMUM_LENGTH+1], sn[MAX_SEQUENCE_NAME_SIZE+1]; 
char infile[MAX_FILE_NAME_SIZE+1]; 
double version=VERSION_NUMBER; 

char alphabet[AA_ALPHABET_SIZE+1] = "ACDEFGHIKLMNPQRSTVWXY";
char sorted_alphabet[AA_ALPHABET_SIZE+1] = "ACDEFGHIKLMNPQRSTVWXY";

bool verbose, too_long, test, headfoot, iset, LR, scrambled_sample; 

double lp[AA_ALPHABET_SIZE]; 
double aa_composition[AA_ALPHABET_SIZE] = { 0.0822, 0.0137, 0.0581, 0.0690, 0.0406, 
                                            0.0719, 0.0230, 0.0589, 0.0586, 0.0943, 
                                            0.0220, 0.0415, 0.0455, 0.0369, 0.0517, 
                                            0.0587, 0.0535, 0.0721, 0.0135, 0.0001, 
                                            0.0343 }; 
                                            /* domains array from fLPS code*/ 

FILE *f4; 


void print_help() 
{ fprintf(stdout, "\nrunny v.%.1lf\n", version);   
  fprintf(stdout,   "==============\n" 
" -h   prints help\n\n"
" -v   verbose\n\n"
" -O   output file prefix; also prefixed are a random number for file uniqueness, and list of options used\n\n"
); 
} /* end of print_help */ 


int file_exist(char *filename) { struct stat sa; return (stat (filename, &sa) == 0); } 


void generate_scrambled_sequence(int seqlen) 
{
char *temp_sequence, temp_char; 

temp_sequence = (char *)calloc(seqlen,sizeof(char)); 

strcpy(temp_sequence,max_sequence);

/* Fisher-Yates shuffle */ 
    for (int i = seqlen - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        temp_char        = temp_sequence[i];
        temp_sequence[i] = temp_sequence[j];
        temp_sequence[j] = temp_char;
    }
    strncpy(randseq, temp_sequence, seqlen);
    randseq[sequence_length] = '\0';

free(temp_sequence); 
} /* end of generate_scrambled_sequence() */ 


int calculate_hpep(const char *seq, int seqlen, bool flag) 
{
int hp=0, aa=-1, previous_aa=-1, current_length=0; 
int i, j, start=0, end=0; 
int *fill; 

fill = (int *)calloc(seqlen,sizeof(int));

if(flag) 
  { residue_hpep = (int *) calloc(alphabet_length,sizeof(int));  
    min_length=MAXIMUM_LENGTH; }

if(!flag)
  {

  i=0; 
  while(i<seqlen-min_length+1)
     { 
     start=i; 
     for(j=i+1;j<seqlen-min_length+1;j++)
        { if(seq[i]==seq[j]) { end=j; i=end; /*fprintf(f4,"loop #1\t%s\t%d\t%d\t%d\n", sequence_name, min_length, start, end);*/ } else { break; } } 
     if(end-start+1>=min_length) 
       { for(j=start;j<=end;j++) 
            { fill[j]=1; } } 
     i++; 
     } /* end of while() */ 

  /*fprintf(f4,"loop #2\t%s\t%d\t%d\t%d\t", sequence_name, min_length, start, end); */
  /*for(j=0;j<seqlen;j++) { fprintf(f4,"%d", fill[j]); } fprintf(f4,"\n"); */
  } /* end of if !flag */ 
else{ /*flag*/ 
    for(i=0;i<seqlen-2;i++)
       { if(seq[i]==seq[i+1] && seq[i+1]==seq[i+2])
           { fill[i]=1; fill[i+1]=1; fill[i+2]=1; } 
       } /* end for i */ 
    } /* end of else flag */ 

for(i=0;i<seqlen;i++)
   { aa=seq[i]; 
   if(fill[i]==1)  
     {
     if(flag) { if(aa==previous_aa) { current_length++; } 
                else { current_length=1; } 
              /*fprintf(f4,"%s\t%d\t%c\t%c\n", sequence_name, current_length, (char) aa, (char) previous_aa); */ 
              } 
     switch(aa) {
          case 65 /*'A'*/:  if(flag) { residue_hpep[0]++; }; hp++; break; 
          case 67 /*'C'*/:  if(flag) { residue_hpep[1]++; }; hp++; break;
          case 68 /*'D'*/:  if(flag) { residue_hpep[2]++; }; hp++; break; 
          case 69 /*'E'*/:  if(flag) { residue_hpep[3]++; }; hp++; break;
          case 70 /*'F'*/:  if(flag) { residue_hpep[4]++; }; hp++; break; 
          case 71 /*'G'*/:  if(flag) { residue_hpep[5]++; }; hp++; break;
          case 72 /*'H'*/:  if(flag) { residue_hpep[6]++; }; hp++; break;
          case 73 /*'I'*/:  if(flag) { residue_hpep[7]++; }; hp++; break;
          case 75 /*'K'*/:  if(flag) { residue_hpep[8]++; }; hp++; break;
          case 76 /*'L'*/:  if(flag) { residue_hpep[9]++; }; hp++; break;
          case 77 /*'M'*/:  if(flag) { residue_hpep[10]++; }; hp++; break;
          case 78 /*'N'*/:  if(flag) { residue_hpep[11]++; }; hp++; break;
          case 80 /*'P'*/:  if(flag) { residue_hpep[12]++; }; hp++; break;
          case 81 /*'Q'*/:  if(flag) { residue_hpep[13]++; }; hp++; break;
          case 82 /*'R'*/:  if(flag) { residue_hpep[14]++; }; hp++; break;
          case 83 /*'S'*/:  if(flag) { residue_hpep[15]++; }; hp++; break;
          case 84 /*'T'*/:  if(flag) { residue_hpep[16]++; }; hp++; break;
          case 86 /*'V'*/:  if(flag) { residue_hpep[17]++; }; hp++; break;
          case 87 /*'W'*/:  if(flag) { residue_hpep[18]++; }; hp++; break;
          case 89 /*'Y'*/:  if(flag) { residue_hpep[20]++; }; hp++; break;
          case 88 /*'X'*/:  if(flag) { residue_hpep[19]++; }; hp++; break; } 
      if(flag) { previous_aa=aa; }  
      } /* end of if fill[i]==1 */ 
   else { /* fill[i]==0 */
        if(flag) { if(current_length<min_length && current_length>0) 
                   { min_length=current_length; current_length=0; } 
                 /*fprintf(f4,"%s\t%d\t%d\t%c\t%c\n", sequence_name, current_length, min_length, (char) aa, (char) previous_aa); */
                 } 
        } /* end of if fill[i]==0 */ 
   } /* end for i */ 

        if(flag) { if(current_length<min_length && current_length>0) 
                   { min_length=current_length; current_length=0; } 
                 /*fprintf(f4,"%s\t%d\t%d\t%c\t%c\n", sequence_name, current_length, min_length, (char) aa, (char) previous_aa); */
                 }    

free(fill); 

return hp; 
} /* end of calculate_hpep() */ 


void convert_aa_sequence()
{
int a,aa; 
char *ps, *pcs; 

ps=sequence; pcs=cleanseq; 
memset(whole_rescounts, 0, alphabet_length * sizeof(int)); 

for(a=0; ps<sequence+sequence_length; ps++)
    { aa=toupper(*ps); 
      switch(aa) {
          case 65 /*'A'*/:  *(pcs+a) = aa; whole_rescounts[0]++; a++; break; 
          case 67 /*'C'*/:  *(pcs+a) = aa; whole_rescounts[1]++; a++; break;
          case 68 /*'D'*/:  *(pcs+a) = aa; whole_rescounts[2]++; a++; break; 
          case 69 /*'E'*/:  *(pcs+a) = aa; whole_rescounts[3]++; a++; break;
          case 70 /*'F'*/:  *(pcs+a) = aa; whole_rescounts[4]++; a++; break; 
          case 71 /*'G'*/:  *(pcs+a) = aa; whole_rescounts[5]++; a++; break;
          case 72 /*'H'*/:  *(pcs+a) = aa; whole_rescounts[6]++; a++; break; 
          case 73 /*'I'*/:  *(pcs+a) = aa; whole_rescounts[7]++; a++; break;
          case 75 /*'K'*/:  *(pcs+a) = aa; whole_rescounts[8]++; a++; break; 
          case 76 /*'L'*/:  *(pcs+a) = aa; whole_rescounts[9]++; a++; break;
          case 77 /*'M'*/:  *(pcs+a) = aa; whole_rescounts[10]++; a++; break; 
          case 78 /*'N'*/:  *(pcs+a) = aa; a++; whole_rescounts[11]++; break;
          case 80 /*'P'*/:  *(pcs+a) = aa; a++; whole_rescounts[12]++; break; 
          case 81 /*'Q'*/:  *(pcs+a) = aa; a++; whole_rescounts[13]++; break;
          case 82 /*'R'*/:  *(pcs+a) = aa; a++; whole_rescounts[14]++; break; 
          case 83 /*'S'*/:  *(pcs+a) = aa; a++; whole_rescounts[15]++; break;
          case 84 /*'T'*/:  *(pcs+a) = aa; a++; whole_rescounts[16]++; break; 
          case 86 /*'V'*/:  *(pcs+a) = aa; a++; whole_rescounts[17]++; break;
          case 87 /*'W'*/:  *(pcs+a) = aa; a++; whole_rescounts[18]++; break; 
          case 89 /*'Y'*/:  *(pcs+a) = aa; a++; whole_rescounts[20]++; break;
          case 88 /*'X'*/:  *(pcs+a) = aa; a++; whole_rescounts[19]++; break; } } 
sequence_length=a; 
} /* end of convert_aa_sequence() */ 


void read_sequence_line(char *sb)
{
int test_length, current_length; 

test_length=strlen(sb);
current_length=strlen(sequence); 
if(test_length+current_length > MAXIMUM_LENGTH)
  { strncat(sequence,sb,MAXIMUM_LENGTH-current_length-1); 
    if(verbose && !too_long) { fprintf(stderr, 
"the length of the sequence %s is too great, only the first %d residues will be analyzed\n", 
sequence_name, MAXIMUM_LENGTH); } 
    too_long=1; 
  } /* end of if test_length+strlen(sequence) > MAXIMUM_LENGTH */ 
else{ strcat(sequence,sb); }
} /* end of read_sequence_line() */ 


void process_sequence()
{
int i,j,k,u; 
int hpep=0, original=0; 
double rn, scrambled_zscore, mean, std, sum, sum2, pvalue; 
double rand_hpep[SCRAMBLED_SAMPLE]; 
char t; 
char hpep_list[500], each_hpep[100]; 

/* make sorted alphabet & max_sequence */ 
max_sequence[0] = '\0'; 
memset(whole_sorted_rescounts, 0, alphabet_length * sizeof(int)); 
memset(rand_hpep, 0.0, SCRAMBLED_SAMPLE * sizeof(int)); 
strcpy(sorted_alphabet,alphabet); 

for(i=0;i<AA_ALPHABET_SIZE;i++)
    { whole_sorted_rescounts[i] = whole_rescounts[i]; } 
            for (i = 1; i < AA_ALPHABET_SIZE; i++) {
                for (j = i; j > 0 && whole_sorted_rescounts[j - 1] > whole_sorted_rescounts[j]; j--) {
                    t = sorted_alphabet[j - 1];
                    sorted_alphabet[j - 1] = sorted_alphabet[j];
                    sorted_alphabet[j] = t;
                    u = whole_sorted_rescounts[j - 1];
                    whole_sorted_rescounts[j - 1] = whole_sorted_rescounts[j];
                    whole_sorted_rescounts[j] = u;  }   } 

sorted_alphabet[AA_ALPHABET_SIZE]='\0'; 
alphabet[AA_ALPHABET_SIZE]='\0'; 

max_sequence_len=0;
for (i=AA_ALPHABET_SIZE-1; i>=0; i--) 
    { if (whole_sorted_rescounts[i]!=0) 
         { for (j = 0; j < whole_sorted_rescounts[i]; j++) 
               { if(max_sequence_len<MAXIMUM_LENGTH -1)
                   { strncat(max_sequence, sorted_alphabet+i,1);
                     max_sequence_len++; }
                 else { fprintf(stderr, 
                            "Warning: max_sequence length exceeded MAXIMUM_LENGTH. Truncating.\n");
                               break; } } /* end for j */ 
                } /* end of whole_sorted_rescounts[i]!=0 */ 
                if (max_sequence_len >= MAXIMUM_LENGTH -1) break;
            } /* end for i */ 
max_sequence[max_sequence_len] = '\0';

/* calculate homopeptide content, pass flag to keep track of residue-specific homopeptide counts */ 
original=1; 
hpep = calculate_hpep(cleanseq, max_sequence_len, original); 

/* scrambled sequences */ 
scrambled_zscore = 9999.11; /* dummy value */ 

              /* SCRAMBLED_SAMPLE */  
              for(k=0,sum=0.0;k<SCRAMBLED_SAMPLE;k++)
                 { generate_scrambled_sequence(max_sequence_len); 
                   rand_hpep[k] = calculate_hpep(randseq, max_sequence_len, 0); 
                   sum += (double) rand_hpep[k]; 
                 } 
 
                  mean = sum / (double) SCRAMBLED_SAMPLE;                 
                  for(k=0,sum2=0.0;k<SCRAMBLED_SAMPLE;k++) 
                     { sum2 += (((double) rand_hpep[k]-mean)*((double) rand_hpep[k]-mean)); } 
                  std=sqrt(sum2/SCRAMBLED_SAMPLE); 
                  scrambled_zscore = (hpep-mean)/std; 

/* output homopeptide content */ 
sprintf(hpep_list,"|"); 
sprintf(each_hpep,""); 
for(i=0;i<alphabet_length;i++)
   { if(residue_hpep[i]>=3) 
       { sprintf(each_hpep,"%c=%d|", alphabet[i], residue_hpep[i]); strcat(hpep_list, each_hpep); }
   } /* end for i */ 

if(erfc(fabs(scrambled_zscore)/sqrt(2.0))>1.0) { pvalue = 1.0; } 
else { pvalue = erfc(fabs(scrambled_zscore)/sqrt(2.0)); } 
if(hpep>0) { 
fprintf(f4, "%s\tRUNNY\t%d\t%d\t%lf\t%d\t%.2le,%.2lf,(%.2lf+/-%.2lf)\t%s\n", 
    sequence_name, max_sequence_len, hpep, ((double) hpep)/((double) max_sequence_len), min_length, 
    pvalue, scrambled_zscore, mean, std, hpep_list); }

} /* end of process_sequence() */ 


void analyse(FILE *fi)
{
bool past_first_name_line=0; 

while(fgets(line,MAXIMUM_LENGTH-1,fi)!=NULL)
{ 
if(!strncmp(line,">",1))
  {
  if(past_first_name_line)
    { 
    sequence_length = strlen(sequence); 
    convert_aa_sequence(); 
    sequence_length = strlen(cleanseq); 
    if(verbose) { fprintf(stderr, "analyzing %s  #%d sequence...\n", sequence_name, number_of_sequences); }

    process_sequence(); 

    sequence[0]='\0'; sequence_name[0]='\0'; cleanseq[0]='\0'; 
    } /* end of if past_first_name_line */ 
  else { past_first_name_line=1; }
  sscanf(line,"%s",sn); strcpy(sequence_name,sn+1); 
  number_of_sequences++; too_long=0; 
  }/*** end of if !strncmp(line,">",1) ***/ 
else { read_sequence_line(line); } 
} /* end of while fgets(line,MAXIMUM_LENGTH-1,f1)!=NULL */ 
if(feof(fi))
  { /* handle last sequence */ 
  sequence_length = strlen(sequence); 
  convert_aa_sequence(); 
  sequence_length = strlen(cleanseq); 
  if(verbose) { fprintf(stderr, "analyzing %s  #%d (last) sequence...\n", sequence_name, number_of_sequences); } 
 
  process_sequence(); 
    
  if(verbose) { fprintf(stderr, "Finished analysis of database %s.\n", infile); } 
  } else { fprintf(stderr, "error in reading database %s, exiting ...\n", infile); exit(1); }
} /* end of analyse_aa() */ 


int main(int argc, char **argv)
  { 
  int i,j; 
  FILE *f1,*f2,*f3; 
  bool NS=0; 
  
  /* for options */ 
  int c, errflg=0, run_id=0, oop=0, opened=0; 
  extern char *optarg;
  extern int optind, optopt; 
  char prefix[MAX_FILE_NAME_SIZE+1]; 
  char outfile[MAX_FILE_NAME_SIZE+1]; 
  char tempopt[50]="", optlist[50]=""; 
  char startrealtime[30]; 
  time_t t; 
  size_t filesize; 

  /* RNG is initialized from the system time */ 
  strcpy(startrealtime,__TIME__); 
  srand(time(&t)); 

/*  *  *  *  PROCESS COMMAND-LINE OPTIONS  *  *  *  */ 
while((c = getopt(argc, argv, "hvli:O:y")) != -1) {
     switch(c) {
     case 'h': print_help(); exit(0);   
     case 'v': verbose=1;     sprintf(tempopt,"%c",optopt); strcat(optlist, tempopt); break;
     /* case 'l': LR=1;     sprintf(tempopt,"%c",optopt); strcat(optlist, tempopt); break;           
     case 'i': iset=1; sscanf(optarg,"%d", &iterations); sprintf(tempopt,"i=%d",optopt); strcat(optlist, tempopt); break; */ 
     case 'O': oop=1; strcpy(prefix,optarg);   sprintf(tempopt,"%c",optopt); strcat(optlist, tempopt); break; 
     /* errors */           
     case 'y': test=1; errflg++; break; 
     case ':': fprintf(stderr, "option -%c requires a value\n", optopt); errflg++; break;
     case '?': fprintf(stderr, "unrecognized option: -%c\n", optopt); errflg++;
} /* end of switch */ 
} /* end of while() getopt */ 
if (errflg) { print_help(); exit(1); } 

/*  *  *  *  OPEN SEQUENCE FILE  *  *  *  */ 
if((f1 = fopen(argv[optind],"r"))==NULL) 
  { fprintf(stderr, "There is no sequence file. Please supply one in FASTA format.\n");
    exit(1); } 
else { strcpy(infile, argv[optind]); } 

/*  *  *  *  HANDLE OUTPUT FILE STREAM(S)  *  *  *  */ 
if(strlen(optlist)==0) { strcpy(optlist,"default"); } 
if(oop) { opened=0;  
while(!opened) { run_id = (int) floor( rand()/100000); 
                 sprintf(outfile,"%s.%d.%s.runny.txt", prefix, run_id, optlist); 
                 if(!file_exist(outfile)) { f4=fopen(outfile, "w"); opened=1; }  
               } /*end of while !opened*/ 
} /* end if oop */ 
else { f4 = stdout; } 

/*  *  *  *  OUTPUT OPTIONS/PARAMETERS, IF VERBOSE  *  *  *  */ 
if(verbose)
{ fprintf(stderr, "%s %s %s\n", argv[0], __DATE__, startrealtime); 
for(i=0;i<argc;i++) { fprintf(stderr, "%s ", argv[i]); } 
} /* end of output parameters if verbose */ 

/*  *  *  *  ALLOCATE GLOBAL ARRAYS *  *  *  */ 
sequence = (char *)calloc(MAXIMUM_LENGTH,sizeof(char)); 
cleanseq = (char *)calloc(MAXIMUM_LENGTH,sizeof(char)); 
randseq  = (char *)calloc(MAXIMUM_LENGTH,sizeof(char)); 
max_sequence = (char *)calloc(MAXIMUM_LENGTH,sizeof(char)); 
whole_rescounts = (int *)calloc(alphabet_length,sizeof(int)); 
whole_sorted_rescounts = (int *)calloc(alphabet_length,sizeof(int)); 

if(!sequence)
  { fprintf(stderr, "memory allocation error, decrease MAXIMUM_LENGTH in code," 
    " ...exiting...\n"); exit(1); }

/* HEADER */ 
fprintf(f4,"NAME/ID\t(RUNNY)\tSEQLEN\tHPEP\tHPEP_PROPORTION\tMIN_TRACT_LENGTH\tNORMALP,ZSCORE,(MEAN+/-STD)\tHPEP_PER_AATYPE\n"); 
fprintf(f4,"=======\t=======\t======\t====\t===============\t================\t===========================\t===============\n"); 

analyse(f1); 

if(verbose) { fprintf(stderr, "runny has finished analysis of %d sequences in database file %s.\n", 
    number_of_sequences, infile); }

exit(0); 
}/******** end of main() ********/ 


