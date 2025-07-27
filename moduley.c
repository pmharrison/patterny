/******************************************************************************/
/*  moduley.c  */ 
/*  code to calculate compositional modules (CModules) and other possible boundaries for */ 
/*  compositionally-biased regions */  
/****  Copyright 2025. Paul Martin Harrison. ****/ 
/****
 ****  Licensed under the 3-clause BSD license. See LICENSE.txt bundled with the package. 
 ****/ 
/****  
 ****  to compile: 
 ****   gcc -O2 -o moduley moduley.c -lm 
 ****
 ****  to run and get help: 
 ****   ./moduley -h 
 **** 
 ****  This program is part of patterny
 ****
 ****  The latest version of this code and patterny generally is available at:
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

#define MAX_REGIONS 100
#define MAX_REGION_LENGTH 5000
#define MAX_SEQUENCE_NAME_SIZE 250 
#define MAXIMUM_LENGTH 15000      
#define MAX_FILE_NAME_SIZE 1000 

double version = 1.0; 

typedef struct {
    char accession[MAX_SEQUENCE_NAME_SIZE];      // Stores the accession string (field $1)
    long tract_start;         // Stores the tract start integer (field $5)
    long tract_end;           // Stores the tract end integer (field $6)
    long bias_count;          // Stores the bias count integer (field $7)
    long seqlen; 
    double bias_pvalue;       // Stores the bias p-value double (field $8)
    char signature[22];       // Stores the signature string (derived from field $9)
    char subsequence[MAX_REGION_LENGTH];    // Stores the subsequence string (field $16)
    int deselect;             // Flag for CMODULES logic, initialized to 0
    int is_redundant;         // Flag for BOUNDARY_SET logic, initialized to 0
} Record;

// Global variables for storing records for the current accession.
// Using a dynamic array to handle an unknown number of records per accession.
Record *records = NULL;
int num_records = 0;        // Current number of records stored
int records_capacity = 0;   // Current allocated capacity for records

// Command-line option flags
bool verbose_mode = false;
bool help_requested = false;
bool fasta_output = false; 

FILE *f4; 

void process_records_for_accession(); 
void make_cmodules();
void make_boundary_set();
void output_cmodules();
void output_boundary_set();
void print_help(); // Function to display help message


const double OVERLAP_THRESHOLD = 0.5;
const int SMALL_MARGIN = 5.0; 

void print_help() {
    fprintf(stderr, "\nmoduley v.%.1lf\n", version);
    fprintf(stderr,
            "=============\n"
            "Usage: moduley [OPTIONS] <input_file>\n\n"
            "Processes biological sequence data to identify modules and boundaries.\n\n"
            "Options:\n"
            "  -h   Display this help message.\n"
            "  -v   Enable verbose output.\n"
            "  -f   Use FASTA output (CModules only in the file).\n"
"  -O   output file prefix; also prefixed are a random number for file uniqueness, and list of options used.\n\n"
            "Arguments:\n"
            "  <input_file>  The path to the input data file.\n\n"
            "Example:\n"
            "  ./moduley -v my_fLPS_olong_output.txt\n\n");
}


int file_exist(char *filename) { struct stat sa; return (stat (filename, &sa) == 0); } 


int main(int argc, char *argv[]) {
    int opt, errflg=0, oop=0, opened=0, run_id=0;
    char prefix[MAX_FILE_NAME_SIZE+1]; 
    char outfile[MAX_FILE_NAME_SIZE+1]; 
    char input_filename[MAX_FILE_NAME_SIZE+1];
    FILE *input_file = NULL;

  char startrealtime[30]; 
  time_t t; 
  size_t filesize; 

  /* RNG is initialized from the system time */ 
  strcpy(startrealtime,__TIME__); 
  srand(time(&t)); 

    // Parse command-line options using getopt
    while ((opt = getopt(argc, argv, "hvfO:")) != -1) {
        switch (opt) {
            case 'h':
                help_requested = true;
                break;
            case 'v':
                verbose_mode = true;
                break;
            case 'f': fasta_output = true; break; 
            case 'O': oop=1; strcpy(prefix,optarg); break; 
            case ':': fprintf(stderr, "option -%c requires a value\n", optopt); errflg++; break;
            case '?': fprintf(stderr, "unrecognized option: -%c\n", optopt); errflg++;
        }
    }
if (errflg) { print_help(); exit(1); } 

    // After getopt, optind points to the first non-option argument.
    // We expect exactly one non-option argument: the input filename.
    if (optind < argc) {
        strcpy(input_filename,argv[optind]);
    } else {
        // No input file provided
        fprintf(stderr, "Error: Input file not specified.\n");
        print_help();
        return EXIT_FAILURE;
    }

    // If help was requested, print it and exit.
    if (help_requested) {
        print_help();
        return EXIT_SUCCESS;
    }

    // Open the input file
    input_file = fopen(input_filename, "r");
    if (input_file == NULL) {
        perror("Error opening input file");
        return EXIT_FAILURE;
    }

    if (verbose_mode) {
        fprintf(stderr, "Verbose mode enabled.\n");
        fprintf(stderr, "Processing file: %s\n", input_filename);
    }


   // HANDLE OUTPUT FILE STREAM(S) 
   if(oop) { opened=0;  
   while(!opened) { run_id = (int) floor( rand()/100000); 
                 sprintf(outfile,"%s.%d.moduley.txt", prefix, run_id); 
                 if(!file_exist(outfile)) { f4=fopen(outfile, "w"); opened=1; }  
               } /*end of while !opened*/ 
   } /* end if oop */ 
   else { f4 = stdout; } 



// HEADER

if(!fasta_output) { 
fprintf(f4,"NAME/ACCESSION\t(MODULEY)\tCMODULES/BOUNDARY_SET\tSEQLEN\tSTART\tEND\tPVALUE\tSIGNATURE\tSUBSEQUENCE\n"); 
fprintf(f4,"==============\t=========\t=====================\t======\t=====\t===\t======\t=========\t===========\n"); } 


    char line[MAXIMUM_LENGTH]; // Buffer to read each line from the input file
    char previous_accession[MAX_SEQUENCE_NAME_SIZE] = "PREVIOUS_START";
    int line_number = 0;

    // Allocate initial memory for records.
    records_capacity = MAX_REGIONS;
    records = (Record *)malloc(records_capacity * sizeof(Record));
    if (records == NULL) {
        perror("Failed to allocate memory for records");
        fclose(input_file); // Close file before exiting on error
        return EXIT_FAILURE;
    }


    // Read input line by line from the specified file until EOF.
    while (fgets(line, sizeof(line), input_file) != NULL) {
        line_number++; // Increment line number for NR > 1 check

        // Temporary variables to hold parsed data from the current line.
        char current_accession[MAX_SEQUENCE_NAME_SIZE];
        long current_tract_start, current_tract_end, current_bias_count, current_seqlen;
        double current_bias_pvalue;
        char current_signature_raw[22]; // To capture the raw $9 before substring
        char current_subsequence[MAX_REGION_LENGTH];


        int fields_read = sscanf(line, "%s\t%ld\t%*s\t%*s\t%ld\t%ld\t%ld\t%le\t%21s\t%*s\t%*s\t%*s\t%*s\t%*s\t%*s\t%s",
                                 current_accession, &current_seqlen, 
                                 &current_tract_start, &current_tract_end,
                                 &current_bias_count, &current_bias_pvalue,
                                 current_signature_raw,
                                 current_subsequence);

        // Check if all expected fields were successfully read.
        // Adjust this number based on the actual number of fields you expect to parse.
        if (fields_read < 7) {
            if (verbose_mode) {
                fprintf(stderr, "Warning: Could not parse all expected fields from line %d: %s", line_number, line);
            }
            continue; // Skip to the next line if parsing fails
        }

        // This block is executed when a new accession is encountered.
        if (strcmp(current_accession, previous_accession) != 0) {
            // If it's not the very first line (NR > 1 in AWK) and there are records
            // from the previous accession, process them.
            if (line_number > 1 && num_records > 0) {
                if (verbose_mode) {
                    fprintf(stderr, "Processing records for accession: %s (previous)\n", previous_accession);
                }
                process_records_for_accession();
            }
            // Reset for the new accession.
            num_records = 0;
            strcpy(previous_accession, current_accession); // Update previous_accession
            if (verbose_mode) {
                fprintf(stderr, "New accession detected: %s\n", current_accession);
            }
        }

        // Add the current record to the dynamic array.
        // If the array is full, reallocate more memory.
        if (num_records >= records_capacity) {
            records_capacity *= 2; // Double the capacity
            records = (Record *)realloc(records, records_capacity * sizeof(Record));
            if (records == NULL) {
                perror("Failed to reallocate memory for records");
                fclose(input_file); // Close file before exiting on error
                return EXIT_FAILURE;
            }
            if (verbose_mode) {
                fprintf(stderr, "Increased records capacity to %d\n", records_capacity);
            }
        }

        // Populate the fields of the current record in the array.
        strcpy(records[num_records].accession, current_accession);
        records[num_records].seqlen = current_seqlen; 
        records[num_records].tract_start = current_tract_start;
        records[num_records].tract_end = current_tract_end;
        records[num_records].bias_count = current_bias_count;
        records[num_records].bias_pvalue = current_bias_pvalue;

        // This extracts the signature string excluding the first and last characters.
        size_t sig_len = strlen(current_signature_raw);
        if (sig_len > 2) { // Ensure there are enough characters to extract
            strncpy(records[num_records].signature, current_signature_raw + 1, sig_len - 2);
            records[num_records].signature[sig_len - 2] = '\0'; // Null-terminate the string
        } else {
            records[num_records].signature[0] = '\0'; // Handle empty or too short signatures
        }

        strcpy(records[num_records].subsequence, current_subsequence);
        records[num_records].deselect = 0;      // Initialize deselect flag
        records[num_records].is_redundant = 0;  // Initialize is_redundant flag
        num_records++; // Increment the count of records for the current accession
    }


    // After reading all lines, process any remaining records for the last accession.
    if (num_records > 0) {
        if (verbose_mode) {
            fprintf(stderr, "Processing remaining records for accession: %s (final)\n", previous_accession);
        }
        process_records_for_accession();
    }


    fclose(input_file); // Close the input file
    free(records); // Free the dynamically allocated memory before exiting.
    return EXIT_SUCCESS;
}


void process_records_for_accession() {
make_cmodules(); 
if(!fasta_output) { make_boundary_set(); }
output_cmodules();
if(!fasta_output) { output_boundary_set(); } 
} /* end of process_records_for_accession() */ 


void make_cmodules() {
    for (int i = 0; i < num_records; i++) {
        if (!records[i].deselect) { // Only process if not already deselected
            for (int j = i + 1; j < num_records; j++) {
                // Checks if signature[i] contains the first character of signature[j].
                char first_char_sig_j[2] = {records[j].signature[0], '\0'}; // Get first char
                if (strlen(records[j].signature) > 0 && strstr(records[i].signature, first_char_sig_j) != NULL) {

                    // If j ends before i starts
                    if (records[j].tract_start < records[i].tract_start &&
                        records[j].tract_end < records[i].tract_start) {
                        continue;
                    }

                    // If i ends before j starts
                    else if (records[i].tract_start < records[j].tract_start &&
                             records[i].tract_end < records[j].tract_start) {
                        continue;
                    }

                    // If i fully contains j
                    else if (records[i].tract_start <= records[j].tract_start &&
                             records[i].tract_end >= records[j].tract_end) {
                        records[j].deselect = 1;
                    }

                    // If j fully contains i
                    else if (records[j].tract_start < records[i].tract_start &&
                             records[j].tract_end > records[i].tract_end) {
                        records[j].deselect = 1;
                    }

                    // If j start is equal 
                    else if (records[j].tract_start == records[i].tract_start) {
                        records[j].deselect = 1;
                    }

                    // If j end is equal 
                    else if (records[j].tract_end == records[i].tract_end ) {
                        records[j].deselect = 1;
                    }

                    // If j overlaps with the start of i
                    else if (records[j].tract_start < records[i].tract_start &&
                             records[j].tract_end >= records[i].tract_start &&
                             records[j].tract_end <= records[i].tract_end) {
                        // Calculate overlap percentage
                        double overlap = (double)(records[j].tract_end - records[i].tract_start + 1) /
                                         (records[j].tract_end - records[j].tract_start + 1);
                        if (overlap >= OVERLAP_THRESHOLD) {
                            records[j].deselect = 1;
                        }
                    }

                    // If j overlaps with the end of i
                    else if (records[j].tract_start <= records[i].tract_end &&
                             records[j].tract_end > records[i].tract_end &&
                             records[j].tract_start >= records[i].tract_start) {
                        // Calculate overlap percentage
                        double overlap = (double)(records[i].tract_end - records[j].tract_start + 1) /
                                         (records[j].tract_end - records[j].tract_start + 1);
                        if (overlap >= OVERLAP_THRESHOLD) {
                            records[j].deselect = 1;
                        }
                    }
                } // End if signature match
            } // End for j loop
        } // End if !deselect[i]
    } // End for i loop
}


void make_boundary_set() {
    for (int i = 0; i < num_records; i++) {
        if (!records[i].is_redundant) { // Only process if not already marked redundant
            for (int j = i + 1; j < num_records; j++) {
                char first_char_sig_j[2] = {records[j].signature[0], '\0'};
                if (strlen(records[j].signature) > 0 && strstr(records[i].signature, first_char_sig_j) != NULL) {

                    // Check if tract start and end are within SMALL_MARGIN of each other
                    // Using fabs for floating-point absolute difference
                    if (fabs((double)records[j].tract_start - (double) records[i].tract_start) <= SMALL_MARGIN &&
                        fabs((double)records[j].tract_end - (double) records[i].tract_end) <= SMALL_MARGIN) {
                        records[j].is_redundant = 1;  } 
                    else if(records[j].tract_start==records[i].tract_start && records[j].tract_end==records[i].tract_end) { 
                        records[j].is_redundant = 1;  } 
                } // End if signature match
            } // End for j loop
        } // End if !is_redundant[i]
    } // End for i loop
}


void output_cmodules() {
int i; 

if(!fasta_output)
{
for(i=0;i<num_records;i++) 
   {
   if(!records[i].deselect) { 
        fprintf(f4, "%s\tMODULEY\tCMODULES\t%ld\t%ld\t%ld\t%.2le\t{%s}\t%s\n",
               records[i].accession,
               records[i].seqlen,
               records[i].tract_start,
               records[i].tract_end,
               records[i].bias_pvalue,
               records[i].signature,
               records[i].subsequence); } /* if !records[i].deselect */ 
    } /* end for i */ 
} /* end of if(!fasta_output) */ 
else { /* fasta_output */ 
for(i=0;i<num_records;i++) 
   {
   if(!records[i].deselect)
     {  
    /* output ">" line and subsequence */ 
    fprintf(f4, ">%s_{%s}_%ld_%ld_%.2le\n", records[i].accession, records[i].signature, 
        records[i].tract_start, records[i].tract_end, records[i].bias_pvalue);  
    fprintf(f4, "%s\n\n", records[i].subsequence); 
     } /* if !records[i].deselect */ 
   } /* end for i*/ 
} /* end of else fasta_output */ 

} /* end of output_cmodules() */ 


void output_boundary_set() {
int i; 
    for(i=0;i<num_records;i++) 
       {
        if(!records[i].is_redundant) { 
        fprintf(f4, "%s\tMODULEY\tBOUNDARY_SET\t%ld\t%ld\t%ld\t%.2le\t{%s}\t%s\n",
               records[i].accession,
               records[i].seqlen,
               records[i].tract_start,
               records[i].tract_end,
               records[i].bias_pvalue,
               records[i].signature,
               records[i].subsequence); } 
    } /* end for i */ 
} /* end of output_cmodules() */ 


