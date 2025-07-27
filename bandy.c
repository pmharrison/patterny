/******************************************************************************/
/*  bandy.c  */ 
/*  code to calculate compositional banding given output from the fLPS program */  
/****  Copyright 2025. Paul Martin Harrison. ****/ 
/****
 ****  Licensed under the 3-clause BSD license. See LICENSE.txt bundled with the package. 
 ****/ 
/****  
 ****  to compile: 
 ****   gcc -O2 -o bandy bandy.c -lm 
 ****
 ****  to run and get help: 
 ****   ./bandy -h 
 **** 
 ****  This program is part of patterny
 ****
 ****  The latest version of this code and patterny generally is available at:
 ****    https://github.com/pmharrison
 **** 
 **** 
 ****    Harrison, PM. "Patterny: A troupe of decipherment helpers for intrinsic disorder, 
 ****     low complexity and compositional bias in proteins", Biomolecules, ... 
 ****/ 
/*****************************************************************************************/ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <ctype.h> 
#include <sys/types.h>
#include <sys/stat.h>   
#include <stdbool.h> 
#include <unistd.h> 
#include <sys/time.h> 
#include <sys/resource.h> 
#include <time.h> 

#define MAX_REGIONS 1000 
#define MAX_REGION_LENGTH 5000 
#define MAXIMUM_LENGTH 15000 
#define MAX_FILE_NAME_SIZE 500 
#define MAX_SEQUENCE_NAME_SIZE 250 
#define MAX_BANDS 1000 
#define MAX_PERFECT_ARRAY 10000 
#define SAMPLE_SIZE 1000 
#define MAD_THRESHOLD 3.5

double version=1.0; 

void print_help() 
{ fprintf(stderr, "\nbandy v.%.1lf\n", version);   
  fprintf(stderr,   "=============\n" 
" -h   prints help\n\n"
" -v   prints verbose info\n\n" 
" -O   ...FILENAME_PREFIX... allows for output to be sent to a unique file prefixed with FILENAME_PREFIX.\n"
"      The parameters that are specified are included in the generated filename.\n\n");   
} /* end of print_help */ 

int file_exist(char *filename) { struct stat sa; return (stat (filename, &sa) == 0); } 

int compare_ints(const void *a, const void *b) { return (*(int*)a - *(int*)b); } 

double calculate_median(int *arr, int n) 
{
int *temp_arr = (int *)malloc(n * sizeof(int));
for(int i=0;i<n;i++) { temp_arr[i] = arr[i]; }
qsort(temp_arr, n, sizeof(int), compare_ints);
double median;
if(n%2==0) { median = (double)(temp_arr[n/2-1] + temp_arr[n/2])/2.0; } 
else { median = (double)temp_arr[n/2]; }
free(temp_arr);
return median;
} /* end of calculate_median() */ 

int main(int argc, char **argv) {
    int f,g,h,i,j,k,l,n,t,u,s,m,w; 
    FILE *f1, *f4; 
    char input_file[MAX_FILE_NAME_SIZE], line[MAXIMUM_LENGTH], output_line[MAXIMUM_LENGTH]; 
    char band_string[MAXIMUM_LENGTH], band_pair[100], outlier_string[MAXIMUM_LENGTH], outlier_pair[100]; 

    char current_bias_type[10], wb; 
    bool overlap; 
    int band_number, clean_band_number, max_band_size, half_max_band_size, highest_end, offset; 
    int p, rand_p[SAMPLE_SIZE], outlier_count=0; 
    double normalp; 

    /* single variables */ 
    char single_accession[MAX_REGIONS][MAX_SEQUENCE_NAME_SIZE], single_signature[MAX_REGIONS][4]; 
    int single_sequence_length[MAX_REGIONS], single_tract_number[MAX_REGIONS], single_tract_start[MAX_REGIONS]; 
    int single_tract_end[MAX_REGIONS], single_bias_count[MAX_REGIONS]; 
    double single_bias_pvalue[MAX_REGIONS]; 

    /* multiple variables */ 
    char multiple_accession[MAX_REGIONS][MAX_SEQUENCE_NAME_SIZE], multiple_signature[MAX_REGIONS][20]; 
    int multiple_sequence_length[MAX_REGIONS], multiple_tract_number[MAX_REGIONS], multiple_tract_start[MAX_REGIONS]; 
    int multiple_tract_end[MAX_REGIONS], multiple_bias_count[MAX_REGIONS]; 
    double multiple_bias_pvalue[MAX_REGIONS]; 

    /* whole variables */ 
    char whole_accession[MAX_REGIONS][MAX_SEQUENCE_NAME_SIZE], whole_signature[MAX_REGIONS][4]; 
    int whole_sequence_length[MAX_REGIONS], whole_tract_number[MAX_REGIONS], whole_tract_start[MAX_REGIONS]; 
    int whole_tract_end[MAX_REGIONS], whole_bias_count[MAX_REGIONS]; 
    double whole_bias_pvalue[MAX_REGIONS]; 

    /* band variables */ 
    int band_start[MAX_BANDS], band_end[MAX_BANDS]; 
    int temp_band_start[MAX_BANDS], temp_band_end[MAX_BANDS], perfect_band_start[MAX_BANDS]; 
    int perfect_band_end[MAX_BANDS], band_centre[MAX_BANDS], perfect_band_size, total_band_size, half_perfect_band_size; 
    int temp_band_number, random_band_start[MAX_BANDS], random_band_end[MAX_BANDS]; 
    int current_band[MAX_BANDS], current_interval[MAX_BANDS], temp_interval[MAX_BANDS]; 
    int random_boundary[2*MAX_BANDS], band_flag[MAX_BANDS];  
    int interval[MAX_BANDS], interval_flag[MAX_BANDS], abs_deviations[MAX_BANDS]; 
    int original_band_start[MAX_BANDS], original_band_end[MAX_BANDS], original_band_number;

    bool no_sampling, pick_again, duplicate; 
    int rb, interval_number; 
    double sum, sum2, mean, std, zscore, median_interval=0.0; 

  /* for options */ 
  int c, errflg=0, run_id=0, oop=0, opened=0; 
  extern char *optarg;
  extern int optind, optopt; 
  char prefix[MAX_FILE_NAME_SIZE+1]; 
  char outfile[MAX_FILE_NAME_SIZE+1]; 
  char tempopt[50]="", optlist[50]=""; 
  char startrealtime[30]; 
  time_t tt; 
  size_t filesize; 
  bool test, verbose=0; 

  /* RNG is initialized from the system time */ 
  strcpy(startrealtime,__TIME__); 
  srand(time(&tt)); 

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

/* input file */ 
if((f1 = fopen(argv[optind],"r"))==NULL)
  { fprintf(stderr, "There is no sequence file. Please supply one in FASTA format.\n");
    exit(1); }   
else { strcpy(input_file, argv[optind]); } 

/*  *  *  *  HANDLE OUTPUT FILE STREAM(S)  *  *  *  */ 
if(strlen(optlist)==0) { strcpy(optlist,"default"); }
if(oop) { opened=0;  
while(!opened) { run_id = (int) floor( rand()/100000); 
                 sprintf(outfile,"%s.%d.%s.bandy.out.txt", prefix, run_id, optlist); 
                 if(!file_exist(outfile)) { f4=fopen(outfile, "w"); opened=1; }  
               } /*end of while !opened*/ 
} /* end if oop */ 
else { f4 = stdout; } 

/*  *  *  *  OUTPUT OPTIONS/PARAMETERS, IF VERBOSE  *  *  *  */ 
if(verbose) 
{ fprintf(f4, "%s %s %s\n", argv[0], __DATE__, startrealtime); 
for(i=0;i<argc;i++) { fprintf(f4, "%s ", argv[i]); } 
fprintf(f4,"\n"); }  /* end of output parameters if verbose */ 

/** header **/ 
if(verbose)
{ fprintf(f4, "(BANDY)\tNAME\tBIAS\tSEQLEN\tREGION#\tBIASCOUNT\tBIASPVALUE\tDPB\tZSCORE\tNORMALP\tBANDCOUNT\tINTERVAL_OUTLIERS:\t\tBANDS:\t\t\n"); 
  fprintf(f4, "=======\t====\t====\t======\t=======\t=========\t==========\t===\t======\t=======\t=========\t==================\t\t======\t\t\n"); }

    /* read the compositional bias data */ 
    n=0; 
    while(fgets(line,MAXIMUM_LENGTH-1,f1)!=NULL)
         { 
         if(!strncmp(line,"#",1)) { s=m=w=0; }

         else { /* read the bias tract data */ 
              if(strncmp(line,"<",1))
                {
                sscanf(line, "%*s %*s %s", current_bias_type); 
                if(!strcmp(current_bias_type,"SINGLE"))
                  { sscanf(line, "%s %d %*s %d %d %d %d %le %s", 
                        single_accession[s], &single_sequence_length[s], &single_tract_number[s], 
                        &single_tract_start[s], &single_tract_end[s], &single_bias_count[s], 
                        &single_bias_pvalue[s], single_signature[s]); 
                    s++; } /* end of if SINGLE */ 

               else if(!strcmp(current_bias_type,"MULTIPLE")) 
                      { sscanf(line, "%s %d %*s %d %d %d %d %le %s", 
                        multiple_accession[m], &multiple_sequence_length[m], &multiple_tract_number[m], 
                        &multiple_tract_start[m], &multiple_tract_end[m], &multiple_bias_count[m], 
                        &multiple_bias_pvalue[m], multiple_signature[m]); 
                        m++; } /* end of if MULTIPLE */ 

               else { /* WHOLE */ 
                      sscanf(line, "%s %d %*s %d %d %d %d %le %s", 
                        whole_accession[w], &whole_sequence_length[w], &whole_tract_number[w], 
                        &whole_tract_start[w], &whole_tract_end[w], &whole_bias_count[w], 
                        &whole_bias_pvalue[w], whole_signature[w]); 
                      w++; } /* end of WHOLE */ 
                } /* end of if not a summary "<" footer */ 
              else { /* if a summary "<" footer */ 
                   /*** process the tract data for the specific sequence ***/ 

                   for(i=0;i<w;i++)
                      { 

                      memset(band_start, 0,  MAX_REGIONS * sizeof(int)); 
                      memset(band_end, 0,  MAX_REGIONS * sizeof(int)); 
                      memset(original_band_start, 0,  MAX_REGIONS * sizeof(int)); 
                      memset(original_band_end, 0,  MAX_REGIONS * sizeof(int));                      
                      band_number=0; 
                      wb=whole_signature[i][1];  

                        for(j=0;j<s;j++)
                           {
                           if(!strncmp(single_signature[j]+1,whole_signature[i]+1,1))
                             { overlap=0; 
                               for(k=0;k<band_number;k++) 
                                  {
if((single_tract_start[j]<=band_start[k]&&single_tract_end[j]>=band_start[k])||(single_tract_start[j]<=band_end[k]&&single_tract_end[j]>=band_end[k])) 
                                   { overlap=1; break; }   
else if(single_tract_start[j]<=band_start[k]&&single_tract_end[j]>=band_start[k]) 
                                   { overlap=1; break; }
else if(single_tract_start[j]<=band_end[k]&&single_tract_end[j]>=band_end[k]) 
                                   { overlap=1; break; }
                                  } /* end for k */ 

                              /* store the potential band */ 
                              if(overlap==0)
                              { band_start[ band_number ] = single_tract_start[j]; 
                                band_end[ band_number ]   = single_tract_end[j]; 
                                band_number++; }

                           } /* end of if !strncmp(single_signature[j]+1,whole_signature[i]+1,1) */ 
                         } /* end for j */ 


                      for(j=0;j<m;j++)
                         {
                         if(!strncmp(multiple_signature[j]+1,whole_signature[i]+1,1))
                              { /* if wb */ 
                              overlap=0; 
                              for(k=0;k<band_number;k++) 
                                 {
if((multiple_tract_start[j]<=band_start[k]&&multiple_tract_end[j]>=band_start[k])||(multiple_tract_start[j]<=band_end[k]&&multiple_tract_end[j]>=band_end[k])) 
                                   { overlap=1; break; }
else if(multiple_tract_start[j]<=band_start[k]&&multiple_tract_end[j]>=band_start[k]) 
                                   { overlap=1; break; }
else if(multiple_tract_start[j]<=band_end[k]&&multiple_tract_end[j]>=band_end[k]) 
                                   { overlap=1; break; }
                                 } /* end for k */ 
                              /* store the potential band */ 
                              if(overlap==0)
                                { band_start[ band_number ] = multiple_tract_start[j]; 
                                  band_end[ band_number ]   = multiple_tract_end[j]; 
                                  band_number++; } 
                              } /* end of if wb */ 
                          } /* end for j */ 


                        /* check for duplications */ 
                        temp_band_number=0; 
                        memset(temp_band_start, 0,  sizeof temp_band_start); 
                        memset(temp_band_end, 0,  sizeof temp_band_end); 
                        for(l=0;l<band_number;l++)
                           {  duplicate=0; 
                              for(k=0;k<l;k++)
                                 { if(band_start[k]==band_start[l] && band_end[k]==band_end[l])
                                     { duplicate=1; break; } 
                                 } /* end for k */ 
                              if(duplicate==0) 
                                { temp_band_start[temp_band_number] = band_start[l]; 
                                  temp_band_end[temp_band_number] = band_end[l]; 
                                  temp_band_number++; }
                                } /* end for l */ 
                        for(l=0;l<temp_band_number;l++) 
                           { band_start[l] = temp_band_start[l]; 
                             band_end[l]   = temp_band_end[l]; }
                             band_number = temp_band_number; 


                        /* sort bands */ 
                        t=u=0; 
                        for(k=1;k<band_number;k++)
                           { for(l=k;l>0&&band_start[l-1]>band_start[l];l--)
                                { t=band_start[l-1]; band_start[l-1]=band_start[l]; band_start[l]=t; 
                                  u=band_end[l-1];   band_end[l-1]=band_end[l];     band_end[l]=u; }  }

                        /* adjust positions to start at zero and store highest end */ 
                        offset = 999999; highest_end = 0; 
                        for(l=0;l<band_number;l++)
                           { if(band_start[l]<offset) 
                               { offset = band_start[l]; } 
                             if(band_end[l]>highest_end) 
                               { highest_end = band_end[l]; }  }

                        total_band_size=0;
                        for(l=0;l<band_number;l++)
                           { original_band_start[l] = band_start[l]; 
                             original_band_end[l]   = band_end[l]; 
                             band_start[l] = band_start[l] - offset; 
                             band_end[l]   = band_end[l]   - offset; 
                             total_band_size += band_end[l]-band_start[l]+1; }

                        highest_end = highest_end - offset ; 
                        original_band_number = band_number; 


                      if(band_number>1) 
                        {
                        if(band_number>3)
                        {
                        /** check for outlier intervals **/ 
                        for(l=0;l<band_number-1;l++) 
                           { interval[l] = band_start[l+1]-band_end[l]-1; } /* just intervening residues */ 
                        memset(interval_flag, 0,  sizeof interval_flag); 
                        memset(abs_deviations, 0,  sizeof abs_deviations); 
                        interval_number=l; 
                        median_interval = calculate_median(interval, interval_number); 
                        for (l=0;l<interval_number;l++) { abs_deviations[l] = (int)fabs((double)interval[l]-median_interval); }
                        double mad_value = calculate_median(abs_deviations, interval_number); 
                        double outlier_threshold_value = MAD_THRESHOLD * mad_value; 

                        int last_outlier_index=0, first_conformer=-1, last_conformer=0;
                        for(l=0;l<interval_number;l++) 
                           { if(fabs((double)interval[l]-median_interval) > outlier_threshold_value) 
                               { outlier_count++; 
                                 if(outlier_count==1) { interval_flag[l]=2; } 
                                 else { interval_flag[l]=1; last_outlier_index=l; } } 

                                 if(interval_flag[l]==0) 
                                   { if(first_conformer<0) { first_conformer=l; }
                                     else { last_conformer=l; } }
                           } /* end for l<interval_number */ 
                        interval_flag[last_outlier_index]=2; 
                        interval_flag[first_conformer]=3; 
                        interval_flag[last_conformer]=3; 

                        /** find outlier intervals **/ 

                        memset(band_flag, 0,  sizeof band_flag); 

                        for(l=0;l<interval_number;l++) 
                           { if(interval_flag[l]==1||interval_flag[l]==2) { band_flag[l]=1; } } 
                        int start_flag_index=-1, end_flag_index=0; 
                        int bandsum=0, mean_pure_band=0; 
                        for(l=0;l<band_number;l++) 
                           { if(band_flag[l]==0)
                               { if(start_flag_index==-1) { start_flag_index=l; }
                                 else { end_flag_index=l; }
                                 bandsum+=(band_end[l]-band_start[l]+1); } 
                           } /* end for l */             
                        /* mean size of pure bands */ 
                        mean_pure_band = (int) ceil((double) bandsum / (double) band_number); 

                        /* search for runs of flagged bands and merge them */ 

                        memset(current_band, 0,  sizeof current_band); 
                        memset(current_interval, 0,  sizeof current_interval); 
                        memset(temp_band_start, 0,  sizeof temp_band_start); 
                        memset(temp_band_end, 0,  sizeof temp_band_end); 
                        memset(temp_interval, 0,  sizeof temp_interval); 
                        int current_band_number=0, merge_counter=0, current_offset=0, run_started=0; 
                        int current_merged_band_start=0, current_merged_band_size=0, current_excised=0, total_excised=0; 
                        int stop_loop=0; 
                        for(l=0;l<band_number-1;l++)
                           { if(!band_flag[l])
                               { current_band[current_band_number]     = band_end[l]-band_start[l]+1; 
                                 current_interval[current_band_number] = interval[l]; 
                                 current_band_number++; }  } 

                        if(!band_flag[band_number-1]) 
                          { current_band[current_band_number] = band_end[band_number-1]-band_start[band_number-1] +1; 
                            current_interval[current_band_number] = 0; 
                            current_band_number++; }                 

                        int temp_band_number=0; 
                        current_offset=0; 
                        for(l=0;l<current_band_number;l++)
                           { temp_band_start[l] = current_offset; 
                             temp_band_end[l] = current_offset + current_band[l] -1; 
                             current_offset += (current_band[l] + current_interval[l]); }

                        /* reassign bands to the band array */ 
                        for(l=0;l<current_band_number;l++)
                           { band_start[l] = temp_band_start[l]; 
                             band_end[l]   = temp_band_end[l]; } 

                        highest_end = band_end[current_band_number-1]; 
                        band_number = current_band_number; 
                        } /* end of if band_number>3 ... check for outlier intervals */ 


                        /* set band parameters */ 
                        max_band_size = (int) ( (double) highest_end / (double) band_number ); 
                        half_max_band_size = (int) ( (double) highest_end / ((double) band_number *2.0));   

                        /* make band_centre array for internal bands */ 
                        band_centre[0] = (int) ((double) max_band_size / 2.0); 
                        band_centre[band_number-1] = (int) highest_end - (int) ((double) max_band_size / 2.0); 
                        for(j=1;j<band_number-1;j++)
                           { band_centre[j] = half_max_band_size + j * max_band_size; }

                        memset(perfect_band_start, 0,  sizeof perfect_band_start); 
                        memset(perfect_band_end,   0,  sizeof perfect_band_end); 

                        perfect_band_size = floor( (double) (total_band_size) / (double) band_number ); 
                        half_perfect_band_size = (int) ((double) perfect_band_size / 2.0); 

                        /* calculate an array of perfect bandings */ 
                        perfect_band_start[0] = 0; 
                        perfect_band_end[0]   = perfect_band_size-1; 
                        perfect_band_start[band_number-1] = highest_end-perfect_band_size+1; 
                        perfect_band_end[band_number-1]   = highest_end; 

                        if(band_number>2) 
                          { for(k=1;k<band_number-1;k++)
                               { perfect_band_start[k] = band_centre[k]-half_perfect_band_size; 
                                 perfect_band_end[k]   = band_centre[k]+half_perfect_band_size; } 
                          } /* end of if band_number>2 */ 

                        /* calculate the distance to perfect banding */ 
                        p=0; 
                        for(l=0;l<band_number;l++)
                           { p += abs(perfect_band_start[l] - band_start[l]); 
                             p += abs(perfect_band_end[l]   - band_end[l]); } 

                        /* calculate 1000 simulated banding of the same band_number
                            |--> calculate the smallest difference to a perfect banding for each; 
                           calculate the mean and SD of these; */ 

                        no_sampling=0; 

                        if(highest_end<=21) /* => not enough room for assortment */ 
                          { 
                            if(verbose) { fprintf(stderr, 
        "simulated banding not possible for sequence %s bias type %s (too short = %d residues)...skipping...\n", 
                            whole_accession[i], whole_signature[i], highest_end); }
                            no_sampling=1; }

                           if(no_sampling==0)
                           {

                           for(h=0;h<SAMPLE_SIZE;h++) 
                              { /** make set of random boundaries across with minimum 1 and maximum highest_end-1 **/ 
                              l=f=0; 

                              for(g=0;g<band_number*2-2;g++) { random_boundary[g]=0; } 

                              while(l<band_number*2-2)
                                   { rb = (int) floor( (double) highest_end * (double) rand()/ (double) RAND_MAX); 
                                     pick_again=0; 
                                     if(rb==0||rb==1||rb==total_band_size-1||rb==total_band_size) { pick_again=1; }
                                     for(g=0;g<l;g++) { if(rb==random_boundary[g]) { pick_again=1; break; }  } 
                                     if(pick_again==0) { random_boundary[l] = rb; l++; } 
                                   } /* end while l */ 

                           /* sort the boundaries numerically */ 
                           t=0; 
                           for(k=1;k<band_number*2-2;k++)
                           { for(l=k;l>0&&random_boundary[l-1]>random_boundary[l];l--)
                                { t=random_boundary[l-1]; random_boundary[l-1]=random_boundary[l]; random_boundary[l]=t; } } 

                           /* assign boundaries to random band arrays */ 
                           random_band_start[0] = 0; 
                           random_band_end[0]   = random_boundary[0]; 
                           random_band_start[band_number-1] = random_boundary[band_number*2-3]; 
                           random_band_end[band_number-1]   = highest_end; 

                           for(g=1,f=1;g<band_number*2-3;g+=2) 
                              { random_band_start[f] = random_boundary[g]; 
                                random_band_end[f]   = random_boundary[g+1]; 
                                f++; }

                           /* calculate p for each random banding and store in array */ 
                           rand_p[h] = 0; 
                           for(l=0;l<band_number;l++)
                              { rand_p[h] += abs(perfect_band_start[l] - random_band_start[l]); 
                                rand_p[h] += abs(perfect_band_end[l]   - random_band_end[l]); } 
                        
                           } /* end of for h */ 

                           /* mean+/-std, zscore */ 
                           sum=sum2=0.0; 
                           for(f=0;f<h;f++) { sum+= (double) rand_p[f]; } 
                           mean= sum / (double) SAMPLE_SIZE; 
                           for(f=0;f<h;f++) { sum2 += (((double) rand_p[f]-mean)*((double) rand_p[f]-mean)); }
                           std=sqrt(sum2/SAMPLE_SIZE); 
                           zscore = ((double)p-mean)/std; 

                           normalp = erfc(fabs(zscore)/sqrt(2.0)); 
                           if(normalp>1.0) { normalp=1.0; }                          

                           strcpy(output_line,""); 
                           strcpy(band_string,"|"); strcpy(band_pair,""); 
                           strcpy(outlier_string,"|"); strcpy(outlier_pair,""); 

                           /**** MODIFY FOR INTERVAL OUTLIERS AND OUTPUT ORIGINAL BANDS LIST ****/
                           sprintf(output_line,"BANDY\t%s\t%s\t%d\t%d\t%d\t%.2le\t%d\t%.2lf\t%.2le\t%d\tOUTLIERS:\t", 
                              /*input_file,*/ whole_accession[i], whole_signature[i], whole_sequence_length[i], whole_tract_number[i], 
                               whole_bias_count[i], whole_bias_pvalue[i], p, zscore, normalp, original_band_number); 

                        /* make outlier_string */ 
                        if(band_number>3) 
                        { if(outlier_count>0) 
                            { for(l=0;l<interval_number;l++)
                                 { if(interval_flag[l]==1||interval_flag[l]==2) 
                                     { sprintf(outlier_pair, "%d:%d|", l+1, interval[l]); 
                                       strcat(outlier_string, outlier_pair); } } 
                            } /* end of if outlier_count>1 */ 
                        else { strcpy(outlier_string,"N/A"); } 
                        } /* end of if band_number>3 */ 
                        else { strcpy(outlier_string,"N/A"); } 

                        /* make band_string */ 
                        for(l=0;l<original_band_number;l++)
                           { sprintf(band_pair, "%d-%d|", original_band_start[l], original_band_end[l]); 
                             strcat(band_string, band_pair); } 
                        strcat(band_string,"\n");  

                        strcat(output_line,outlier_string); 
                        strcat(output_line,"\tBANDS:\t"); 
                        strcat(output_line,band_string); 
                        fputs(output_line,f4); 
                        
                        } /* if no_sampling==0 */ 

                        else { /* no_sampling==1 */  
                        fprintf(f4, "BANDY\t%s\t%s\t%d\t%d\t%d\t%.2le\tN/A\tN/A\tN/A\tN/A\tOUTLIERS: N/A\tBANDS:\tN/A\t\n",
                            /*input_file,*/ whole_accession[i], whole_signature[i], whole_sequence_length[i], whole_tract_number[i], 
                               whole_bias_count[i], whole_bias_pvalue[i]); }

                        } /* end of if band_number>1 */ 
                      } /* end for i<w */ 
                   s=m=w=0; 
                   } /* end of if footer */ 
              } /* end of read bias tract data */  
         } /* end of while */ 
 
if(verbose) { fprintf(stderr, "Finished analyzing for banding.\n"); } 
} /* end of main() */ 

