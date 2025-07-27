/******************************************************************************/
/*  blocky.c  */ 
/*  code to calculate the degree to which amino acids are segregated by type 
    into blocks in input sequences  */  
/****  Copyright 2025. Paul Martin Harrison. ****/ 
/****
 ****  Licensed under the 3-clause BSD license. See LICENSE.txt bundled with the package. 
 ****/ 
/****  
 ****  to compile: 
 ****   gcc -O2 -fno-common -o blocky blocky.c -lm 
 ****
 ****  to run and get help: 
 ****   ./blocky -h 
 **** 
 ****  This program is part of patterny
 ****
 ****  The latest version of this code is available at:
 ****    https://github.com/pmharrison
 **** 
 ****  Citations: 
 ****    Harrison, PM. "Intrinsically Disordered Compositional Bias in Proteins: Sequence Traits, 
 ****     Region Clustering, and Generation of Hypothetical Functional Associations." 
 ****     Bioinform Biol Insights (2024) 18: 11779322241287485. 
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
#define MAXIMUM_LENGTH 100000
#define VERSION_NUMBER 1.0
#define MAX_SEQUENCE_NAME_SIZE 250 
#define MAX_FILE_NAME_SIZE 500 

#define ITERATIONS 25000 
#define SCRAMBLED_SAMPLE 1000 
#define MAXIMUM_SEQUENCES 250000 

#include "blocky_tables.h"

int alphabet_length=AA_ALPHABET_SIZE, sequence_length=MAXIMUM_LENGTH, number_of_sequences; 

char *max_sequence, *min_blocky; 
int *whole_rescounts, *whole_sorted_rescounts; 

char *sequence, *cleanseq, *randseq, sequence_name[MAX_SEQUENCE_NAME_SIZE+1] = {"initial"}; 
char line[MAXIMUM_LENGTH+1], sn[MAX_SEQUENCE_NAME_SIZE+1]; 
char infile[MAX_FILE_NAME_SIZE+1]; 
double version=VERSION_NUMBER; 

char alphabet[AA_ALPHABET_SIZE+1] = "ACDEFGHIKLMNPQRSTVWXY";
char sorted_alphabet[AA_ALPHABET_SIZE+1] = "ACDEFGHIKLMNPQRSTVWXY"; 

int max_sequence_len, iterations=ITERATIONS;
int s; 
char stored_sequence_name[MAXIMUM_SEQUENCES][MAX_SEQUENCE_NAME_SIZE]; 
char stored_top_biases[MAXIMUM_SEQUENCES][500];
double stored_log_ratio[MAXIMUM_SEQUENCES]; 
double stored_log_length[MAXIMUM_SEQUENCES]; 
double stored_blockiness[MAXIMUM_SEQUENCES]; 
double stored_scrambled_zscore[MAXIMUM_SEQUENCES],   stored_mean[MAXIMUM_SEQUENCES],   stored_std[MAXIMUM_SEQUENCES]; 
double stored_scrambled_zscore_R[MAXIMUM_SEQUENCES], stored_mean_R[MAXIMUM_SEQUENCES], stored_std_R[MAXIMUM_SEQUENCES]; 

int residue_blockiness_s[AA_ALPHABET_SIZE]; 
int residue_blockiness_d[AA_ALPHABET_SIZE]; 
double residue_blockiness_r[AA_ALPHABET_SIZE]; 
double randb[SCRAMBLED_SAMPLE], randr[SCRAMBLED_SAMPLE];

bool verbose, too_long, test, headfoot, iset, LR, scrambled_sample; 

double slope, intercept, pearson_r, p_value; 

double lp[AA_ALPHABET_SIZE]; 
double aa_composition[AA_ALPHABET_SIZE] = { 0.0822, 0.0137, 0.0581, 0.0690, 0.0406, 0.0719, 0.0230, 0.0589, 0.0586, 0.0943, 0.0220, 
                                            0.0415, 0.0455, 0.0369, 0.0517, 0.0587, 0.0535, 0.0721, 0.0135, 0.0001, 0.0343 }; 
                                            /* domains array from fLPS code*/ 

FILE *f4; 


void print_help() 
{ fprintf(stdout, "\nblocky v.%.1lf\n", version);   
  fprintf(stdout,   "==============\n" 
" -h   prints help\n\n"
" -v   verbose\n\n"
" -O   output file prefix; also prefixed are a random number for file uniqueness, and list of options used\n\n"
); 
} /* end of print_help */ 


int file_exist(char *filename) { struct stat sa; return (stat (filename, &sa) == 0); } 


double logbinomial(int n, int k, double expect)
{ return logfactorial[n] - logfactorial[n-k] - logfactorial[k] + ((double)k)*log10(expect) + ((double)(n-k))*log10(1.0-expect); } 


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
} /* end generate_scrambled_sequence() */ 


double calculate_blockiness(const char *seq, const int count[], int sequence_length) 
{
double d = 0.0, s = 0.0; 
int i, current_residue_index = 0; 

memset(residue_blockiness_d, 0, alphabet_length * sizeof(int)); 
memset(residue_blockiness_s, 0, alphabet_length * sizeof(int)); 
memset(residue_blockiness_r, 0.0, alphabet_length * sizeof(double)); 

    for (i = 0; i < sequence_length; i++) {
        int valid = 0;
        int min_s = sequence_length, min_d = sequence_length;
        for (int h = 0; h < AA_ALPHABET_SIZE; h++) {
            if (seq[i] == sorted_alphabet[h] && count[h] > 0) { 
                current_residue_index = h; 
                valid = 1;
                break;
            } 
        } 
        if (valid) { 
            for (int j = 0; j < sequence_length; j++) { 
                if (j != i) {
                    int current_interval = abs(j - i);
                    if (seq[i] != seq[j] && current_interval < min_d) {
                        min_d = current_interval;
                    } else if (seq[i] == seq[j] && current_interval < min_s) {
                        min_s = current_interval;
                    } 
                } /* end of if j != i */ 
            } /* end for j */ 

        } /* end of if valid */ 

        /* keep array of d and s for each residue */ 
        residue_blockiness_d[current_residue_index] += min_d; 
        residue_blockiness_s[current_residue_index] += min_s; 
        d += (double) min_d; 
        s += (double) min_s; 
    } /* end for i */ 

    for(i=0;i<AA_ALPHABET_SIZE;i++)
       { if(count[i]>0){ residue_blockiness_r[i] = ((double) residue_blockiness_d[i])/ ((double) residue_blockiness_s[i]); } ; }

return d / s;
} /* end of calculate_blockiness() */ 


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
"the length of the sequence %s is too great, only the first %d residues will be analyzed\n", sequence_name, MAXIMUM_LENGTH); } 
    too_long=1; 
  } /* end of if test_length+strlen(sequence) > MAXIMUM_LENGTH */ 
else{ strcat(sequence,sb); }
} /* end of read_sequence_line() */ 


void process_sequence()
{
int i,j,k,u, non_zero=0; 
char t, primary_bias; 
double original_blockiness, max_blockiness, min_blockiness, initial_min_blockiness; 
double rn; 
int valid_iterations, first_position, second_position, last_non_zero; 
double ratio, p;
double Bmean, Bstd, Bsum, Bsum2; 
double Rmean, Rstd, Rsum, Rsum2; 
double blockiness, scrambled_zscore, scrambled_zscore_R; 
char top_biases[500], each_bias[100]; 


max_sequence[0] = '\0';
non_zero=0; 

memset(whole_sorted_rescounts, 0, alphabet_length * sizeof(int)); 
strcpy(sorted_alphabet,alphabet); 

    for(i=0;i<21;i++)
       { if(whole_rescounts[i]>0) { non_zero++; last_non_zero=whole_rescounts[i]; } 
         whole_sorted_rescounts[i] = whole_rescounts[i]; }

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

    if(non_zero>1) 
      {

    /* MAKE A MAX_SEQUENCE */ 
    max_sequence_len=0;
    for (i=AA_ALPHABET_SIZE-1; i>=0; i--) {
                if (whole_sorted_rescounts[i]!=0) {
                    for (j = 0; j < whole_sorted_rescounts[i]; j++) {
                        if(max_sequence_len<MAXIMUM_LENGTH -1)
                        { strncat(max_sequence, sorted_alphabet+i,1);
                          max_sequence_len++; }
                        else { fprintf(stderr, 
                                "Warning: max_sequence length exceeded MAXIMUM_LENGTH. Truncating.\n");
                               break; } } /* end for j */ 
                } /* end of whole_sorted_rescounts[i]!=0 */ 
                if (max_sequence_len >= MAXIMUM_LENGTH -1) break;
            } /* end for i */ 
            max_sequence[max_sequence_len] = '\0';

    original_blockiness = calculate_blockiness(cleanseq, whole_sorted_rescounts, max_sequence_len);

    /* calculate the binomial probabilities from whole_sorted_rescounts */ 
    sprintf(top_biases,"|"); 
    sprintf(each_bias,""); 
    for(i=AA_ALPHABET_SIZE-1;i>=0;i--) 
       { if(whole_sorted_rescounts[i]>1 && sorted_alphabet[i]!='X')
           { lp[i] = logbinomial(max_sequence_len, whole_sorted_rescounts[i], aa_composition[i]); 
             sprintf(each_bias, "%c,%d,%.1e,%.3lf|", 
               sorted_alphabet[i], whole_sorted_rescounts[i], pow(10.0,lp[i]), residue_blockiness_r[i]); 
             strcat(top_biases, each_bias); } } 

            sequence_length = max_sequence_len; 
            scrambled_zscore = 9999.11; /* dummy value */ 

              /* SCRAMBLED_SAMPLE */  
              for(k=0,Bsum=0.0,Rsum=0.0;k<SCRAMBLED_SAMPLE;k++)
                 { generate_scrambled_sequence(max_sequence_len); 
                   randb[k] = calculate_blockiness(randseq, whole_sorted_rescounts, sequence_length); 
                   Bsum += randb[k]; } 
 
                  Bmean = Bsum / (double) SCRAMBLED_SAMPLE;  

                  for(k=0,Bsum2=0.0,Rsum2=0.0;k<SCRAMBLED_SAMPLE;k++) 
                     { Bsum2 += (((double) randb[k]-Bmean)*((double) randb[k]-Bmean)); } 

                  Bstd=sqrt(Bsum2/SCRAMBLED_SAMPLE); 
                  scrambled_zscore = (original_blockiness-Bmean)/Bstd; 

                  if(erfc(fabs(scrambled_zscore)/sqrt(2.0))>1.0) { p=1.0; } 
                  else { p = erfc(fabs(scrambled_zscore)/sqrt(2.0)); }

              /* output */ 
                   fprintf(f4, "%s\t%d\t%.2lf\t%.2le,%.2lf,(%.2lf+/-%.2lf)\t%s\n", 
                     sequence_name, max_sequence_len, original_blockiness, 
                      p, scrambled_zscore, Bmean, Bstd, top_biases); 

        } else if(non_zero==1) { fprintf(f4, "%s\t%d\t%.2lf\t%.4lf\thomopeptide\n", 
                                   sequence_name, last_non_zero, 1.0, 1.0); 
                               } /* end else if non_zero == 1 */ 

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
  double p=1.0; 
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
while((c = getopt(argc, argv, "hvO:y")) != -1) {
     switch(c) {
     case 'h': print_help(); exit(0);   
     case 'v': verbose=1;     sprintf(tempopt,"%c",optopt); strcat(optlist, tempopt); break;
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
                 sprintf(outfile,"%s.%d.%s.blocky.txt", prefix, run_id, optlist); 
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
  { fprintf(stderr, "memory allocation error, decrease MAXIMUM_LENGTH and UNIT_STORE in code," 
    " ...exiting...\n"); exit(1); }

/*  *  *  *  ANALYZE SEQUENCES *  *  *  */ 

/* header for output */ 
fprintf(f4, "NAME/ACCESSION\tSEQLEN\tB(BLOCKINESS)\tNORMALP,ZSCORE,(MEAN_B +/- STD)\tBLOCKINESS_FOR_AATYPES:|AA,COUNT,BINOMIALP,B(BLOCKINESS)|...\n"); 
fprintf(f4, "==============\t======\t=============\t===============================\t============================================================\n\n"); 

analyse(f1); 

if(verbose) 
{ fprintf(stderr, "blocky has finished analysis of %d sequences in database file %s.\n", number_of_sequences, infile); } 

exit(0); 
}/******** end of main() ********/ 

