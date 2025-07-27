/******************************************************************************/
/*  pickbanding.c  */ 
/*  program to pick out compositional bandings according to various criteria */    
/****  Copyright 2025. Paul Martin Harrison. ****/ 
/****
 ****  Licensed under the 3-clause BSD license. See LICENSE.txt bundled with the package. 
 ****/ 
/****  
 ****  to compile: 
 ****   gcc -O2 -o pickbanding pickbanding.c -lm 
 ****  
 ****  to run and get help: 
 ****   ./pickbanding -h 
 **** 
 ****  This program is part of patterny
 ****
 ****  The latest version of this code and patterny generally is available at:
 ****    https://github.com/pmharrison
 **** 
 ****  Citations: 
 **** 
 ****    Harrison, PM. "Patterny: A troupe of decipherment helpers for intrinsic disorder, 
 ****     low complexity and compositional bias in proteins", Biomolecules, ... 
 ****/ 
/*****************************************************************************************/ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h> // For DBL_MAX and -DBL_MAX (though not strictly used for initialization, good for reference)

#define MAXIMUM_LENGTH 10000 
#define MAX_FILE_NAME_SIZE 500 
#define MAX_SEQUENCE_NAME_SIZE 250 

// Define a structure to hold the data from each line
typedef struct {
    char original_line[MAXIMUM_LENGTH]; // To store the entire line
    char field2[MAX_SEQUENCE_NAME_SIZE];
    char field3[20];
    double field9;
    int field11;
    int is_valid; // Flag to indicate if parsing was successful and numerical fields are valid
} LineData;

// Function to parse a line and extract relevant fields
// Returns 1 on success, 0 on failure (e.g., N/A values or not enough fields)
int parseLine(char *line, LineData *data) {
    char *token;
    char *rest = line;
    int field_count = 0;
    char temp_line[MAXIMUM_LENGTH]; // Use a temporary buffer for tokenization
    strncpy(temp_line, line, sizeof(temp_line) - 1);
    temp_line[sizeof(temp_line) - 1] = '\0'; // Ensure null-termination

    // Store the original line
    strncpy(data->original_line, line, sizeof(data->original_line) - 1);
    data->original_line[sizeof(data->original_line) - 1] = '\0';

    data->is_valid = 1; // Assume valid until proven otherwise

    // Tokenize the line by tab delimiter
    while ((token = strsep(&rest, "\t")) != NULL) {
        field_count++;
        // Extract field 2
        if (field_count == 2) {
            strncpy(data->field2, token, sizeof(data->field2) - 1);
            data->field2[sizeof(data->field2) - 1] = '\0';
        }
        // Extract field 3
        else if (field_count == 3) {
            strncpy(data->field3, token, sizeof(data->field3) - 1);
            data->field3[sizeof(data->field3) - 1] = '\0';
        }
        // Extract field 9 (double)
        else if (field_count == 9) {
            if (strcmp(token, "N/A") == 0) {
                data->is_valid = 0; // Mark as invalid
                return 0;
            }
            data->field9 = atof(token);
        }
        // Extract field 11 (integer)
        else if (field_count == 11) {
            if (strcmp(token, "N/A") == 0) {
                data->is_valid = 0; // Mark as invalid
                return 0;
            }
            data->field11 = atoi(token);
        }
    }
    // Ensure all required fields were found (at least 11 fields)
    if (field_count < 11) {
        data->is_valid = 0;
        return 0;
    }
    return 1;
}

// Function to initialize a LineData struct to an invalid state
void initializeLineData(LineData *data) {
    data->is_valid = 0;
    data->field9 = 0.0; // Default value, will be overwritten by first valid line
    data->field11 = 0; // Default value
    data->original_line[0] = '\0';
    data->field2[0] = '\0';
    data->field3[0] = '\0';
}

// Function to copy line data
void copyLineData(LineData *dest, const LineData *src) {
    *dest = *src; // Simple struct copy
}


int main(int argc, char **argv) {
    FILE *file;
    char line[MAXIMUM_LENGTH], input_file[MAX_FILE_NAME_SIZE];

    // Pointers to store the best lines for each criterion for the current group
    LineData *best_abs9_then_high11 = NULL; // Criterion (i)
    LineData *best_high9 = NULL;            // Criterion (ii)
    LineData *best_low9 = NULL;             // Criterion (iii)

    // Variables to hold the group identifiers
    char current_group_field2[MAXIMUM_LENGTH] = "";
    char current_group_field3[MAXIMUM_LENGTH] = "";
    int first_line_processed = 0;


    // input file 
    if((file = fopen(argv[1],"r"))==NULL)
      { fprintf(stderr, "There is no input file for pickbanding. Please supply one.\n");
        exit(1); }   
    else { strcpy(input_file, argv[1]); } 


    // Allocate memory for the best lines
    best_abs9_then_high11 = (LineData *)malloc(sizeof(LineData));
    best_high9 = (LineData *)malloc(sizeof(LineData));
    best_low9 = (LineData *)malloc(sizeof(LineData));

    if (!best_abs9_then_high11 || !best_high9 || !best_low9) {
        perror("Memory allocation failed");
        if (best_abs9_then_high11) free(best_abs9_then_high11);
        if (best_high9) free(best_high9);
        if (best_low9) free(best_low9);
        fclose(file);
        return 1;
    }

    // Initialize the best lines to an invalid state
    initializeLineData(best_abs9_then_high11);
    initializeLineData(best_high9);
    initializeLineData(best_low9);


printf("OUTPUT_CRITERION\t(BANDY)\tNAME/ACCESSION\t{PRIMARY_BIAS}\tSEQLEN\tTRACT_NUMBER\tPRIMARY_BIAS_COUNT\tPRIMARY_BIAS_PVALUE\tDPB\tZSCORE\tNORMALP\tBAND_COUNT\tOUTLIER:_INTERVAL_LIST\tBAND:_LIST\n"); 
printf("================\t=======\t==============\t==============\t======\t============\t==================\t===================\t===\t======\t=======\t==========\t======================\t==========\n"); 


    // Process the file line by line
    while (fgets(line, sizeof(line), file) != NULL) {
        // Remove trailing newline character if present
        line[strcspn(line, "\n")] = 0; // Remove newline
        LineData current_parsed_line;
        if (!parseLine(line, &current_parsed_line)) {
            continue; // Skip lines with N/A or invalid data
        }

        // Initialize group with the first valid line encountered
        if (!first_line_processed) {
            strncpy(current_group_field2, current_parsed_line.field2, sizeof(current_group_field2) - 1);
            strncpy(current_group_field3, current_parsed_line.field3, sizeof(current_group_field3) - 1);
            current_group_field2[sizeof(current_group_field2) - 1] = '\0';
            current_group_field3[sizeof(current_group_field3) - 1] = '\0';

            copyLineData(best_abs9_then_high11, &current_parsed_line);
            copyLineData(best_high9, &current_parsed_line);
            copyLineData(best_low9, &current_parsed_line);
            first_line_processed = 1;
            continue; // Move to the next line
        }

        // Check if the current line belongs to the same group
        if (strcmp(current_parsed_line.field2, current_group_field2) == 0 &&
            strcmp(current_parsed_line.field3, current_group_field3) == 0) {

            // Update for Criterion (i): highest 11th field, then highest absolute 9th field
            if (current_parsed_line.field11 > best_abs9_then_high11->field11) {
                copyLineData(best_abs9_then_high11, &current_parsed_line);
            } else if (current_parsed_line.field11 == best_abs9_then_high11->field11) {
                // If 11th field values are equal, compare absolute 9th field
                if (fabs(current_parsed_line.field9) > fabs(best_abs9_then_high11->field9)) {
                    copyLineData(best_abs9_then_high11, &current_parsed_line);
                }
            }

            // Update for Criterion (ii): highest (most positive) 9th field
            if (current_parsed_line.field9 > best_high9->field9) {
                copyLineData(best_high9, &current_parsed_line);
            }

            // Update for Criterion (iii): lowest (most negative) 9th field
            if (current_parsed_line.field9 < best_low9->field9) {
                copyLineData(best_low9, &current_parsed_line);
            }

        } else {
            // New group detected, output results for the previous group
            /* printf("--- Group: %s %s ---\n", current_group_field2, current_group_field3); */ 
            
        if (best_abs9_then_high11->is_valid) {
            printf("HIGHEST_BAND_NUMBER:\t%s\n", best_abs9_then_high11->original_line);
        } else {
            printf("HIGHEST_BAND_NUMBER:\t%s\t%s\tN/A\n", current_group_field2, current_group_field3);
        }
        if (best_high9->is_valid) {
            printf("HIGHEST_ZSCORE:\t%s\n", best_high9->original_line);
        } else {
            printf("HIGHEST_ZSCORE:\t%s\t%s\tN/A\n", current_group_field2, current_group_field3);
        }
        if (best_low9->is_valid) {
            printf("LOWEST_ZSCORE:\t%s\n", best_low9->original_line);
        } else {
            printf("HIGHEST_ZSCORE:\t%s\t%s\tN/A\n", current_group_field2, current_group_field3);
        }
            printf("\n"); // Add a newline for separation

            // Start a new group with the current line
            strncpy(current_group_field2, current_parsed_line.field2, sizeof(current_group_field2) - 1);
            strncpy(current_group_field3, current_parsed_line.field3, sizeof(current_group_field3) - 1);
            current_group_field2[sizeof(current_group_field2) - 1] = '\0';
            current_group_field3[sizeof(current_group_field3) - 1] = '\0';

            copyLineData(best_abs9_then_high11, &current_parsed_line);
            copyLineData(best_high9, &current_parsed_line);
            copyLineData(best_low9, &current_parsed_line);
        }
    }

    // Output results for the last group
    if (first_line_processed) { // Check if any group was processed

        /* printf("--- Group: %s %s ---\n", current_group_field2, current_group_field3); */ 

        if (best_abs9_then_high11->is_valid) {
            printf("HIGHEST_BAND_NUMBER:\t%s\n", best_abs9_then_high11->original_line);
        } else {
            printf("HIGHEST_BAND_NUMBER:\t%s\t%s\tN/A\n", current_group_field2, current_group_field3);
        }
        if (best_high9->is_valid) {
            printf("HIGHEST_ZSCORE:\t%s\n", best_high9->original_line);
        } else {
            printf("HIGHEST_ZSCORE:\t%s\t%s\tN/A\n", current_group_field2, current_group_field3);
        }
        if (best_low9->is_valid) {
            printf("LOWEST_ZSCORE:\t%s\n", best_low9->original_line);
        } else {
            printf("HIGHEST_ZSCORE:\t%s\t%s\tN/A\n", current_group_field2, current_group_field3);
        }
    }

    // Clean up
    free(best_abs9_then_high11);
    free(best_high9);
    free(best_low9);
    fclose(file);

    return 0;
}
