/******************************************************************************/
/*  repeaty.c  */ 
/*  code to calculate repetitiveness of protein tracts */ 
/****  Copyright 2025. Paul Martin Harrison. ****/ 
/****
 ****  Licensed under the 3-clause BSD license. See LICENSE.txt bundled with the package. 
 ****/ 
/****  
 ****  to compile: 
 ****   gcc -O2 -o repeaty repeaty.c -lm 
 ****
 ****  to run and get help: 
 ****   ./repeaty -h 
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
#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdbool.h> // For bool type
#include <errno.h> // For strerror

// Constants
#define MAX_SEQ_LEN 5000
#define MAX_INTERVAL_LEN 100 /* This is the absolute maximum possible internal length for arrays */
#define DEFAULT_MAX_OUTPUT_INTERVAL_LEN 100 /* Default max length for output */
#define ALPHABET_SIZE 256    /* Keep 256 for ASCII indexing for mapping table */
#define INITIAL_POSITIONS_CAPACITY 4
#define N_SCRAMBLES_MAX 1000 // Absolute maximum for N_SCRAMBLES
#define SIGNIFICANCE_THRESHOLD 0.05
#define MAX_FILE_NAME_SIZE 500
#define MAX_SEQUENCE_NAME_SIZE 250
#define NUM_STANDARD_AA 20 /* This is descriptive, actual used for alphabet is VALID_AA_COUNT */
#define MAX_LINE_LENGTH 5000
#define MAX_AA_SYMBOLS 30

/* Global flags for command-line options */
int verbose = 0;
FILE *f4; /* Main output file stream (stdout by default or specified by -O) */

// These variables are now fixed to their default values as command-line options are removed.
const int min_count_threshold = 3; /* Fixed minimum count for an interval to be outputted */
const int max_output_interval_len = DEFAULT_MAX_OUTPUT_INTERVAL_LEN; /* Fixed maximum interval length to output */
int num_scrambles_to_run = 500; /* Number of scrambles, can be changed by -s */
const int sort_by_pvalue = 1; /* Fixed to sort by p-value */
const int long_output_format = 0; /* Fixed to short output format */


/* Define the valid amino acid alphabet including 'X' */
const char* VALID_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWXY";
/* Calculate the actual number of valid amino acids dynamically */
int VALID_AA_COUNT; /* Will be initialized in main */

/* Global array to map amino acid chars to compact indices (0-VALID_AA_COUNT-1) */
int aa_char_to_mapped_index[ALPHABET_SIZE];

// Mathematical helper for Normal CDF (restored from original)
double erf_approx(double x) {
    double a1 = 0.254829592;
    double a2 = -0.284496736;
    double a3 = 1.421413741;
    double a4 = -1.453152027;
    double a5 = 1.061405429;
    double p = 0.3275911;

    int sign = (x < 0) ? -1 : 1;
    x = fabs(x);

    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);
    return sign * y;
}

double normal_cdf(double z) {
    return 0.5 * (1.0 + erf_approx(z / sqrt(2.0)));
}

/* Structure to store a single start-end position pair */
typedef struct {
    int start_pos;
    int end_pos;
} PositionPair;

/* Structure to store all positions for a specific (len, a, y) interval */
typedef struct {
    PositionPair* positions; /* Dynamic array of all occurrences */
    size_t count;            /* Number of occurrences */
    size_t capacity;         /* Allocated capacity */
} AllIntervalPositions;

/* Structure to store scrambled entropy values for a specific interval type */
typedef struct {
    double* plogp_values; /* Dynamic array of p*log(p) values from scrambles */
    size_t count;         /* Number of stored values */
    size_t capacity;         /* Allocated capacity */
    double sum_plogp;     /* Sum for mean calculation */
    double sum_sq_plogp;  /* Sum of squares for variance calculation */
} ScrambledEntropyStats;


/* Structure to hold aggregated data for each interval type (for sorting and summary) */
typedef struct {
    char seq_id[MAX_SEQUENCE_NAME_SIZE+ 1]; /* To store the ID of the sequence */
    char a_char; /* Storing original char, not mapped index */
    char y_char; /* Storing original char, not mapped index */
    int length;
    int count; /* This count is the aggregated total for this specific a-y-len type */
    double plogp_value; /* Observed p*log(p) (now -p*log2(p)) */
    double mean_scrambled_plogp;
    double std_dev_scrambled_plogp;
    double z_score;
    double p_value; /* One-tailed p-value (assuming we are interested in more negative values) */
    /* New fields for SameLength_-plog2p and its components */
    double samelength_plog2p; /* -p*log2(p) calculated from intervals of the same length */
    double samelength_p;      /* p calculated from intervals of the same length */
} IntervalData;

/* Structure to hold the three main entropy values for a sequence */
typedef struct {
    double total_entropy;
    double same_residue_entropy;
    double diff_residue_entropy;
} Entropies;


/* Comparison function for qsort - now always sorts by P-value */
int compareIntervalData(const void* a, const void* b) {
    const IntervalData* ia = (const IntervalData*)a;
    const IntervalData* ib = (const IntervalData*)b;

    /* Primary sort: increasing P-value */
    if (ia->p_value < ib->p_value) return -1;
    if (ia->p_value > ib->p_value) return 1;

    /* Secondary sort: decreasing plogp_value if P-values are equal */
    if (ia->plogp_value > ib->plogp_value) return -1;
    if (ia->plogp_value < ib->plogp_value) return 1;
    
    return 0;
}

/* Comparison function for qsort to sort by count (decreasing) */
int compareIntervalDataByCount(const void* a, const void* b) {
    const IntervalData* ia = (const IntervalData*)a;
    const IntervalData* ib = (const IntervalData*)b;

    if (ia->count > ib->count) return -1; /* Higher count first */
    if (ia->count < ib->count) return 1;

    /* Secondary sort: by P-value (increasing) if counts are equal */
    if (ia->p_value < ib->p_value) return -1;
    if (ia->p_value > ib->p_value) return 1;

    return 0;
}

/* Function to print help message - UPDATED */
void print_help(char *prog_name) {
    printf("Usage: %s [OPTIONS] <fasta_file>\n", prog_name);
    printf("Analyze intervals between amino acids in protein sequences.\n\n");
    printf("Output Format:\n");
    printf("  Repeaty outputs a short summary for each sequence:\n");
    printf("    - The top 10 significant intervals sorted by P-value.\n");
    printf("    - The top 10 significant intervals sorted by Count (P-value as tie-breaker).\n");
    printf("    - A summary line.\n\n");
    printf("Options:\n");
    printf("  -h          Display this help message and exit.\n");
    printf("  -v          Enable verbose output (prints processing messages to stderr).\n");
    printf("  -O <prefix> Specify an output file prefix. Output will be <prefix>.<run_id>.<options>.repeaty.txt\n");
    printf("  -s <num>    Number of scrambles to run (default: %d).\n", num_scrambles_to_run);
}


int file_exist(char *filename) { struct stat sa; return (stat (filename, &sa) == 0); }

/* Function to initialize the character to index mapping */
void initialize_aa_mapping() {
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        aa_char_to_mapped_index[i] = -1; /* Initialize all to -1 (invalid) */
    }
    for (int i = 0; i < VALID_AA_COUNT; i++) {
        aa_char_to_mapped_index[(unsigned char)VALID_AMINO_ACIDS[i]] = i;
    }
}

/* Function to get the mapped index for a given amino acid character */
int get_aa_mapped_index(char c) {
    return aa_char_to_mapped_index[(unsigned char)toupper((unsigned char)c)];
}

/* Function to check if a character is a valid amino acid */
bool is_valid_amino_acid(char c) {
    return get_aa_mapped_index(c) != -1;
}

/* Function to count occurrences of each amino acid */
/* This function counts the frequency of each uppercase amino acid character */
/* in the given sequence, storing counts in the aa_counts array indexed by ASCII value. */
void count_amino_acids(const char* sequence, int* aa_counts) {
    /* Initialize all counts to zero */
    for (int i = 0; i < ALPHABET_SIZE; i++) { /* aa_counts still uses full ASCII range for simplicity */
        aa_counts[i] = 0;
    }
    int seq_len = strlen(sequence);
    /* Iterate through the sequence and count uppercase amino acids that are valid */
    for (int i = 0; i < seq_len; i++) {
        if (is_valid_amino_acid(sequence[i])) {
            aa_counts[(unsigned char)toupper((unsigned char)sequence[i])]++;
        }
    }
}

/* Function to check if an interval should be considered */
/* This function determines if an interval, defined by characters 'a' and 'y', */
/* should be included in the analysis. The criteria are that both 'a' and 'y' */
/* must be valid amino acids (from VALID_AMINO_ACIDS) and their counts in the sequence (from aa_counts) */
/* must be at least 3. */
int should_consider_interval(char a, char y, const int* aa_counts) {
    return is_valid_amino_acid(a) && is_valid_amino_acid(y) &&
           aa_counts[(unsigned char)toupper((unsigned char)a)] >= min_count_threshold && aa_counts[(unsigned char)toupper((unsigned char)y)] >= min_count_threshold;
}

/* Initialize an AllIntervalPositions struct */
void init_all_interval_positions(AllIntervalPositions* aip) {
    aip->positions = NULL;
    aip->count = 0;
    aip->capacity = 0;
}

/* Add a position pair to AllIntervalPositions */
/* Dynamically resizes the positions array if needed. */
int add_position(AllIntervalPositions* aip, int start, int end) {
    if (aip->count == aip->capacity) {
        size_t new_capacity = (aip->capacity == 0) ? INITIAL_POSITIONS_CAPACITY : aip->capacity * 2;
        PositionPair* new_positions = (PositionPair*)realloc(aip->positions, new_capacity * sizeof(PositionPair));
        if (new_positions == NULL) {
            return 0; /* Allocation failed */
        }
        aip->positions = new_positions;
        aip->capacity = new_capacity;
    }
    aip->positions[aip->count].start_pos = start;
    aip->positions[aip->count].end_pos = end;
    aip->count++;
    return 1; /* Success */
}

/* Allocate the 3D interval_positions_map on the heap */
/* Uses VALID_AA_COUNT for the 2nd and 3rd dimensions for memory efficiency */
AllIntervalPositions*** allocate_interval_positions_map(int max_len, int valid_aa_count) {
    AllIntervalPositions*** arr;
    arr = (AllIntervalPositions***)calloc((max_len + 1), sizeof(AllIntervalPositions**));
    if (arr == NULL) return NULL;

    for (int i = 0; i <= max_len; i++) {
        arr[i] = (AllIntervalPositions**)calloc(valid_aa_count, sizeof(AllIntervalPositions*));
        if (arr[i] == NULL) {
            for (int k = 0; k < i; k++) free(arr[k]);
            free(arr);
            return NULL;
        }
        for (int j = 0; j < valid_aa_count; j++) {
            arr[i][j] = (AllIntervalPositions*)calloc(valid_aa_count, sizeof(AllIntervalPositions));
            if (arr[i][j] == NULL) {
                for (int l = 0; l < j; l++) free(arr[i][l]);
                for (int k = 0; k < i; k++) {
                    for (int l = 0; l < valid_aa_count; l++) free(arr[k][l]);
                    free(arr[i]);
                    free(arr);
                    return NULL;
                }
            }
            /* Initialize each AllIntervalPositions struct */
            for (int k = 0; k < valid_aa_count; k++) {
                init_all_interval_positions(&(arr[i][j][k]));
            }
        }
    }
    return arr;
}

/* Function to free the 3D interval_positions_map */
void free_interval_positions_map(AllIntervalPositions*** arr, int max_len, int valid_aa_count) {
    if (arr == NULL) return;
    for (int i = 0; i <= max_len; i++) {
        if (arr[i] != NULL) {
            for (int j = 0; j < valid_aa_count; j++) {
                if (arr[i][j] != NULL) {
                    for (int k = 0; k < valid_aa_count; k++) {
                        free(arr[i][j][k].positions); /* Free inner dynamic array */
                    }
                    free(arr[i][j]);
                }
            }
            free(arr[i]);
        }
    }
    free(arr);
}

/* Initialize ScrambledEntropyStats struct */
void init_scrambled_entropy_stats(ScrambledEntropyStats* ses) {
    ses->plogp_values = (double*)malloc(num_scrambles_to_run * sizeof(double)); // Pre-allocate to num_scrambles_to_run
    if (ses->plogp_values == NULL) {
        // Handle allocation failure, though for N_SCRAMBLES=1000 this should be rare
        fprintf(stderr, "Repeaty: Error: Failed to pre-allocate plogp_values. Please use a computer with more available RAM or reduce the number of scrambled sequences using the -s option.\n");
        ses->capacity = 0;
    } else {
        ses->capacity = num_scrambles_to_run;
    }
    ses->count = 0;
    ses->sum_plogp = 0.0;
    ses->sum_sq_plogp = 0.0;
}

/* Add a plogp value to ScrambledEntropyStats */
int add_scrambled_plogp(ScrambledEntropyStats* ses, double plogp) {
    // If initial pre-allocation failed or a larger capacity is needed (shouldn't for num_scrambles_to_run if pre-allocated enough)
    if (ses->count == ses->capacity) {
        // This case should ideally not happen if num_scrambles_to_run is correctly used for initial capacity.
        // If it does, it indicates a logic error or an unexpectedly high number of scrambles.
        // For robustness, we'll realloc, but it's a warning sign.
        size_t new_capacity = (ses->capacity == 0) ? num_scrambles_to_run : ses->capacity * 2; // Use num_scrambles_to_run as base
        double* new_values = (double*)realloc(ses->plogp_values, new_capacity * sizeof(double));
        if (new_values == NULL) {
            fprintf(stderr, "Repeaty: Error: Failed to reallocate plogp_values for ScrambledEntropyStats.\n");
            return 0; /* Allocation failed */
        }
        ses->plogp_values = new_values;
        ses->capacity = new_capacity;
    }
    ses->plogp_values[ses->count] = plogp;
    ses->sum_plogp += plogp;
    ses->sum_sq_plogp += plogp * plogp;
    ses->count++;
    return 1; /* Success */
}

/* Allocate the 3D scrambled_entropy_map on the heap */
/* Uses VALID_AA_COUNT for the 2nd and 3rd dimensions for memory efficiency */
ScrambledEntropyStats*** allocate_scrambled_entropy_map(int max_len, int valid_aa_count) {
    ScrambledEntropyStats*** arr;
    arr = (ScrambledEntropyStats***)calloc((max_len + 1), sizeof(ScrambledEntropyStats**));
    if (arr == NULL) return NULL;

    for (int i = 0; i <= max_len; i++) {
        arr[i] = (ScrambledEntropyStats**)calloc(valid_aa_count, sizeof(ScrambledEntropyStats*));
        if (arr[i] == NULL) {
            for (int k = 0; k < i; k++) free(arr[k]);
            free(arr);
            return NULL;
        }
        for (int j = 0; j < valid_aa_count; j++) {
            arr[i][j] = (ScrambledEntropyStats*)calloc(valid_aa_count, sizeof(ScrambledEntropyStats));
            if (arr[i][j] == NULL) {
                for (int l = 0; l < j; l++) free(arr[i][l]);
                for (int k = 0; k < i; k++) {
                    for (int l = 0; l < valid_aa_count; l++) free(arr[k][l]);
                    free(arr[i]);
                    free(arr);
                    return NULL;
                }
            }
            for (int k = 0; k < valid_aa_count; k++) {
                init_scrambled_entropy_stats(&(arr[i][j][k]));
                if (arr[i][j][k].plogp_values == NULL) { // Check if init_scrambled_entropy_stats failed
                    // Propagate error cleanup if inner allocation failed
                    for (int x = 0; x <= i; ++x) {
                        if (arr[x] != NULL) {
                            for (int y = 0; y < valid_aa_count; ++y) {
                                if (arr[x][y] != NULL) {
                                    for (int z = 0; z < valid_aa_count; ++z) {
                                        if (x < i || (x == i && y < j) || (x == i && y == j && z < k)) {
                                            free(arr[x][y][z].plogp_values);
                                        }
                                    }
                                    free(arr[x][y]);
                                }
                            }
                            free(arr[x]);
                        }
                    }
                    free(arr);
                    return NULL;
                }
            }
        }
    }
    return arr;
}

/* Free the 3D scrambled_entropy_map */
void free_scrambled_entropy_map(ScrambledEntropyStats*** arr, int max_len, int valid_aa_count) {
    if (arr == NULL) return;
    for (int i = 0; i <= max_len; i++) {
        if (arr[i] != NULL) {
            for (int j = 0; j < valid_aa_count; j++) {
                if (arr[i][j] != NULL) {
                    for (int k = 0; k < valid_aa_count; k++) {
                        free(arr[i][j][k].plogp_values);
                    }
                    free(arr[i][j]);
                }
            }
            free(arr[i]);
        }
    }
    free(arr);
}

/* Generate scrambled sequence (Fisher-Yates shuffle) based on amino acid counts */
void generate_scrambled_sequence(const int* original_aa_counts, int sequence_length, char* scrambled_sequence) {
    char temp_pool[MAX_SEQ_LEN + 1];
    int current_pos = 0;

    /* Populate temp_pool with characters based on original_aa_counts for valid amino acids */
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        /* Only add if it's a valid amino acid and has occurrences */
        if (original_aa_counts[i] > 0 && is_valid_amino_acid((char)i)) {
            for (int j = 0; j < original_aa_counts[i]; j++) {
                if (current_pos < sequence_length) { // Ensure we don't exceed the original sequence length
                    temp_pool[current_pos++] = (char)i;
                } else {
                     // This warning is less critical now as we cap current_pos at sequence_length
                     // if the sum of counts exceeds it (which it shouldn't if aa_counts is from the sequence).
                     break;
                }
            }
        }
    }
    temp_pool[current_pos] = '\0'; /* Null-terminate the pool */

    /* Adjust sequence_length to the actual number of valid amino acids added to the pool */
    /* This ensures we only shuffle the valid AAs present in the original sequence */
    int actual_shuffle_length = current_pos;

    /* Fisher-Yates shuffle */
    for (int i = actual_shuffle_length - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        char temp_char = temp_pool[i];
        temp_pool[i] = temp_pool[j];
        temp_pool[j] = temp_char;
    }
    strncpy(scrambled_sequence, temp_pool, actual_shuffle_length);
    scrambled_sequence[actual_shuffle_length] = '\0';

    // If the original sequence had non-AA characters, pad with 'X' or similar if needed
    // For now, we assume the scrambled_sequence will be used up to its actual_shuffle_length
    // and the rest of the buffer is implicitly handled by analysis functions.
    // If the original sequence had non-AA characters that contributed to seq_len,
    // this scrambled_sequence might be shorter. The analysis functions need to handle this.
    // Given the `is_valid_amino_acid` filter in `read_fasta_sequence`, this should match.
}

/* Function to allocate the 3D interval_counts array on the heap */
/* Uses VALID_AA_COUNT for the 2nd and 3rd dimensions for memory efficiency */
int*** allocate_interval_counts(int max_len, int valid_aa_count) {
    int*** arr;
    arr = (int***)calloc((max_len + 1), sizeof(int**));
    if (arr == NULL) return NULL;

    for (int i = 0; i <= max_len; i++) {
        arr[i] = (int**)calloc(valid_aa_count, sizeof(int*));
        if (arr[i] == NULL) {
            for (int k = 0; k < i; k++) free(arr[k]);
            free(arr);
            return NULL;
        }
        for (int j = 0; j < valid_aa_count; j++) {
            arr[i][j] = (int*)calloc(valid_aa_count, sizeof(int));
            if (arr[i][j] == NULL) {
                for (int l = 0; l < j; l++) free(arr[i][l]);
                for (int k = 0; k < i; k++) {
                    for (int l = 0; l < valid_aa_count; l++) free(arr[k][l]);
                    free(arr[i]);
                    free(arr);
                    return NULL;
                }
            }
            /* Initialize the allocated memory to zero */
            memset(arr[i][j], 0, valid_aa_count * sizeof(int));
        }
    }
    return arr;
}

/* Function to free the 3D interval_counts array */
void free_interval_counts(int*** arr, int max_len, int valid_aa_count) {
    if (arr == NULL) return;
    for (int i = 0; i <= max_len; i++) {
        if (arr[i] != NULL) {
            for (int j = 0; j < valid_aa_count; j++) {
                free(arr[i][j]);
            }
            free(arr[i]);
        }
    }
    free(arr);
}


/* Function to perform interval analysis on a given sequence */
/* This function is designed to be called for both original and scrambled sequences. */
/* It returns the total interval count. If scramble_stats_map is provided, it populates */
/* it with plogp values from the current sequence. It also returns the three summary */
/* entropy values via the calculated_entropies_out pointer. */
int analyze_sequence_intervals(const char* current_sequence, int seq_len,
                               const int* aa_counts, int max_interval_len_for_seq,
                               ScrambledEntropyStats*** scramble_stats_map,
                               Entropies* calculated_entropies_out,
                               int*** current_interval_counts_in) { // Added parameter for pre-allocated counts

    int interval_count = 0;
    // Use memset to efficiently clear the pre-allocated array instead of re-allocating
    for (int i = 0; i <= max_interval_len_for_seq; i++) { // Clear only up to current max_interval_len_for_seq
        for (int j = 0; j < VALID_AA_COUNT; j++) {
            memset(current_interval_counts_in[i][j], 0, VALID_AA_COUNT * sizeof(int));
        }
    }


    /* Find and count intervals for the current sequence (original or scrambled) */
    for (int i = 0; i < seq_len; i++) {
        for (int j = i + 1; j < seq_len; j++) {
            int interval_len = j - i - 1; // Length of sequence *between* i and j
            if (interval_len >= 0 && interval_len <= max_interval_len_for_seq) {
                // Only consider intervals if both start and end AAs meet the minimum count threshold
                if (should_consider_interval(current_sequence[i], current_sequence[j], aa_counts)) {
                    int mapped_a = get_aa_mapped_index(current_sequence[i]);
                    int mapped_y = get_aa_mapped_index(current_sequence[j]);
                    if (mapped_a != -1 && mapped_y != -1) { // Should always be true due to should_consider_interval
                         current_interval_counts_in[interval_len][mapped_a][mapped_y]++;
                         interval_count++;
                    }
                }
            }
        }
    }

    /* Initialize summary entropies for this sequence */
    double current_total_entropy = 0.0;
    double current_same_residue_entropy = 0.0;
    double current_diff_residue_entropy = 0.0;

    if (interval_count > 0) {
        for (int len = 0; len <= max_interval_len_for_seq; len++) {
            for (int a_idx = 0; a_idx < VALID_AA_COUNT; a_idx++) {
                for (int y_idx = 0; y_idx < VALID_AA_COUNT; y_idx++) {
                    if (current_interval_counts_in[len][a_idx][y_idx] > 0) {
                        double p = (double)current_interval_counts_in[len][a_idx][y_idx] / interval_count;
                        /* Negated entropy calculation */
                        double plogp = -p * log2(p);

                        /* If scramble_stats_map is provided, add to the aggregated scrambled stats for individual intervals */
                        if (scramble_stats_map != NULL) {
                            if (!add_scrambled_plogp(&scramble_stats_map[len][a_idx][y_idx], plogp)) {
                                fprintf(stderr, "Repeaty: Error: Failed to add scrambled plogp value for len %d, a %c, y %c.\n",
                                        len, VALID_AMINO_ACIDS[a_idx], VALID_AMINO_ACIDS[y_idx]);
                                return 0; /* Indicate failure */
                            }
                        }

                        /* Accumulate for current sequence's summary entropies */
                        current_total_entropy += plogp;
                        if (a_idx == y_idx) { // Comparing mapped indices
                            current_same_residue_entropy += plogp;
                        } else {
                            current_diff_residue_entropy += plogp;
                        }
                    } else if (scramble_stats_map != NULL) {
                        // If count is 0 for this interval type in a scramble, add 0.0 to maintain array size
                        if (!add_scrambled_plogp(&scramble_stats_map[len][a_idx][y_idx], 0.0)) {
                            fprintf(stderr, "Repeaty: Error: Failed to add zero scrambled plogp value for len %d, a %c, y %c.\n",
                                    len, VALID_AMINO_ACIDS[a_idx], VALID_AMINO_ACIDS[y_idx]);
                            return 0; /* Indicate failure */
                        }
                    }
                }
            }
        }
    } else { // If interval_count is 0 for this sequence (e.g., very short sequence or no valid AAs)
        if (scramble_stats_map != NULL) {
            // Fill all plogp_values with 0.0 for this scramble to maintain array size
            for (int len = 0; len <= max_interval_len_for_seq; len++) {
                for (int a_idx = 0; a_idx < VALID_AA_COUNT; a_idx++) {
                    for (int y_idx = 0; y_idx < VALID_AA_COUNT; y_idx++) {
                        if (!add_scrambled_plogp(&scramble_stats_map[len][a_idx][y_idx], 0.0)) {
                            fprintf(stderr, "Repeaty: Error: Failed to add zero scrambled plogp value (total 0 intervals).\n");
                            return 0; /* Indicate failure */
                        }
                    }
                }
            }
        }
    }

    /* Pass back calculated summary entropies if pointer is valid */
    if (calculated_entropies_out != NULL) {
        calculated_entropies_out->total_entropy = current_total_entropy;
        calculated_entropies_out->same_residue_entropy = current_same_residue_entropy;
        calculated_entropies_out->diff_residue_entropy = current_diff_residue_entropy;
    }

    return interval_count;
}

/* Helper function to print header for short interval output */
void print_short_interval_header(FILE* f_out) {
    fprintf(f_out, "SeqID\tInterval\tPositions\tLength\tCount\t-plog2p\tScrambleMean\tScrambleStdDev\tZ-score\tP-value\n");
}

/* Helper function to print interval details for the short output format (top 10 lists) */
void print_short_interval_line(FILE* f_out, const IntervalData* interval, AllIntervalPositions*** original_interval_positions_map) {
    char a_char_val = interval->a_char;
    char y_char_val = interval->y_char;
    int len_val = interval->length;

    int mapped_a = get_aa_mapped_index(a_char_val);
    int mapped_y = get_aa_mapped_index(y_char_val);

    AllIntervalPositions* current_aip = NULL;
    // Use MAX_INTERVAL_LEN for bounds check here as original_interval_positions_map is allocated with it
    if (mapped_a != -1 && mapped_y != -1 && len_val >= 0 && len_val <= MAX_INTERVAL_LEN) {
        current_aip = &original_interval_positions_map[len_val][mapped_a][mapped_y];
    }

    char pos_str[MAX_SEQ_LEN * 10]; // Max possible positions can lead to a very long string
    pos_str[0] = '\0';
    if (current_aip != NULL && current_aip->count > 0) {
        size_t current_str_len = 0;
        for (size_t p_idx = 0; p_idx < current_aip->count; p_idx++) {
            char temp_pos[50]; // Buffer for one position pair (e.g., "12345-67890,")
            int needed_chars = snprintf(temp_pos, sizeof(temp_pos), "%d-%d%s",
                                        current_aip->positions[p_idx].start_pos + 1,
                                        current_aip->positions[p_idx].end_pos + 1,
                                        (p_idx < current_aip->count - 1 ? "," : ""));

            if (needed_chars < 0) {
                fprintf(stderr, "Repeaty: Warning: Error during position string formatting for interval %c...%c length %d.\n", a_char_val, y_char_val, len_val);
                if (pos_str[0] == '\0') strncpy(pos_str, "Error", sizeof(pos_str) - 1);
                break;
            }
            // Check if appending will exceed buffer size
            if (current_str_len + needed_chars + 1 > sizeof(pos_str)) {
                if (current_str_len + 3 + 1 <= sizeof(pos_str)) { // Make space for "..."
                    strcat(pos_str, "...");
                }
                break;
            }
            strcat(pos_str, temp_pos);
            current_str_len += needed_chars;
        }
    } else {
        strncpy(pos_str, "N/A", sizeof(pos_str) - 1);
        pos_str[sizeof(pos_str) - 1] = '\0';
    }

    fprintf(f_out, "%s\t%c...%c\t%s\t%d\t%d\t%.2e\t%.2e\t%.2e\t%.2f\t%.2e\n",
           interval->seq_id,
           interval->a_char, interval->y_char,
           pos_str, // Added back
           interval->length,
           interval->count,
           interval->plogp_value,
           interval->mean_scrambled_plogp,
           interval->std_dev_scrambled_plogp,
           interval->z_score,
           interval->p_value);
}


/* Helper function to process a single sequence and output results */
void process_sequence_data(const char* header, const char* sequence,
                           FILE* f4_out,
                           int verbose_flag, int min_count, int max_len_output) {

    int seq_len = strlen(sequence);
    // max_interval_len_for_seq_current is now (seq_len / 2)
    // This variable is no longer used for allocation sizes of the main 3D arrays.
    // It's kept here for the loop bounds to ensure we don't process intervals longer than half the sequence,
    // and for the cap at MAX_INTERVAL_LEN.
    int max_interval_len_for_seq_current = (seq_len > 0) ? (seq_len / 2) : 0;
    // Cap at MAX_INTERVAL_LEN to prevent excessive processing if seq_len is very large
    if (max_interval_len_for_seq_current > MAX_INTERVAL_LEN) {
        max_interval_len_for_seq_current = MAX_INTERVAL_LEN;
    }
    if (max_interval_len_for_seq_current < 0) {
        max_interval_len_for_seq_current = 0;
    }

    // Amino acid counts for the original sequence (used for scrambling and filtering)
    int aa_counts[ALPHABET_SIZE];
    count_amino_acids(sequence, aa_counts);

    // Allocate memory for the data structures specific to this sequence's processing
    // These now use MAX_INTERVAL_LEN for their first dimension
    AllIntervalPositions*** original_interval_positions_map = allocate_interval_positions_map(MAX_INTERVAL_LEN, VALID_AA_COUNT);
    if (original_interval_positions_map == NULL) {
        fprintf(stderr, "Repeaty: Error: Failed to allocate memory for original interval positions map for sequence %s. Skipping.\n", header);
        return;
    }

    // Allocate current_interval_counts once here for reuse in analyze_sequence_intervals
    int*** current_interval_counts = allocate_interval_counts(MAX_INTERVAL_LEN, VALID_AA_COUNT);
    if (current_interval_counts == NULL) {
        fprintf(stderr, "Repeaty: Error: Failed to allocate memory for current interval counts for sequence %s. Skipping.\n", header);
        free_interval_positions_map(original_interval_positions_map, MAX_INTERVAL_LEN, VALID_AA_COUNT);
        return;
    }

    // --- 1. Analyze the Original Sequence ---
    Entropies original_sequence_entropies;
    // analyze_sequence_intervals still uses max_interval_len_for_seq_current for loop bounds
    int original_interval_total_count = analyze_sequence_intervals(sequence, seq_len, aa_counts, max_interval_len_for_seq_current, NULL, &original_sequence_entropies, current_interval_counts);

    double observed_total_entropy = original_sequence_entropies.total_entropy;
    double observed_same_residue_entropy = original_sequence_entropies.same_residue_entropy;
    double observed_diff_residue_entropy = original_sequence_entropies.diff_residue_entropy;

    // Populate original_interval_positions_map
    if (original_interval_total_count > 0) { // Only if there are intervals to process
        for (int i_seq = 0; i_seq < seq_len; i_seq++) {
            for (int j_seq = i_seq + 1; j_seq < seq_len; j_seq++) {
                int interval_len = j_seq - i_seq - 1;
                if (interval_len >= 0 && interval_len <= max_interval_len_for_seq_current) { // Use current for loop bounds
                    if (should_consider_interval(sequence[i_seq], sequence[j_seq], aa_counts)) {
                        int mapped_a = get_aa_mapped_index(sequence[i_seq]);
                        int mapped_y = get_aa_mapped_index(sequence[j_seq]);
                        if (mapped_a != -1 && mapped_y != -1) {
                            if (!add_position(&original_interval_positions_map[interval_len][mapped_a][mapped_y], i_seq, j_seq)) {
                                fprintf(stderr, "Repeaty: Error: Failed to add position to original interval map for %s. Skipping.\n", header);
                                free_interval_positions_map(original_interval_positions_map, MAX_INTERVAL_LEN, VALID_AA_COUNT);
                                free_interval_counts(current_interval_counts, MAX_INTERVAL_LEN, VALID_AA_COUNT);
                                return;
                            }
                        }
                    }
                }
            }
        }
    }


    // Calculate total intervals at each length for same-length entropy calculation
    // These arrays are now sized by MAX_INTERVAL_LEN
    double* sum_samelength_plog2p_by_length = (double*)calloc(MAX_INTERVAL_LEN + 1, sizeof(double));
    if (sum_samelength_plog2p_by_length == NULL) {
        fprintf(stderr, "Repeaty: Error: Failed to allocate memory for sum_samelength_plog2p_by_length for sequence %s. Skipping.\n", header);
        free_interval_positions_map(original_interval_positions_map, MAX_INTERVAL_LEN, VALID_AA_COUNT);
        free_interval_counts(current_interval_counts, MAX_INTERVAL_LEN, VALID_AA_COUNT);
        return;
    }

    int* total_intervals_at_length = (int*)calloc(MAX_INTERVAL_LEN + 1, sizeof(int));
    if (total_intervals_at_length == NULL) {
        fprintf(stderr, "Repeaty: Error: Failed to allocate memory for total_intervals_at_length for sequence %s. Skipping.\n", header);
        free(sum_samelength_plog2p_by_length);
        free_interval_positions_map(original_interval_positions_map, MAX_INTERVAL_LEN, VALID_AA_COUNT);
        free_interval_counts(current_interval_counts, MAX_INTERVAL_LEN, VALID_AA_COUNT);
        return;
    }


    for (int len_calc = 0; len_calc <= max_interval_len_for_seq_current; len_calc++) { // Loop up to current max
        for (int a_idx_calc = 0; a_idx_calc < VALID_AA_COUNT; a_idx_calc++) {
            for (int y_idx_calc = 0; y_idx_calc < VALID_AA_COUNT; y_idx_calc++) {
                total_intervals_at_length[len_calc] += original_interval_positions_map[len_calc][a_idx_calc][y_idx_calc].count;
            }
        }
    }

    for (int len_calc = 0; len_calc <= max_interval_len_for_seq_current; len_calc++) { // Loop up to current max
        for (int a_idx_calc = 0; a_idx_calc < VALID_AA_COUNT; a_idx_calc++) {
            for (int y_idx_calc = 0; y_idx_calc < VALID_AA_COUNT; y_idx_calc++) {
                int original_count_for_type = original_interval_positions_map[len_calc][a_idx_calc][y_idx_calc].count;
                if (original_count_for_type > 0) {
                    double p_same_length_temp = 0.0;
                    if (total_intervals_at_length[len_calc] > 0) {
                        p_same_length_temp = (double)original_count_for_type / total_intervals_at_length[len_calc];
                        if (p_same_length_temp > 0) {
                            sum_samelength_plog2p_by_length[len_calc] += -p_same_length_temp * log2(p_same_length_temp);
                        }
                    }
                }
            }
        }
    }
    
    // --- 2. Analyze Scrambled Sequences ---
    // This now uses MAX_INTERVAL_LEN for allocation
    ScrambledEntropyStats*** scrambled_entropy_map = allocate_scrambled_entropy_map(MAX_INTERVAL_LEN, VALID_AA_COUNT);
    if (scrambled_entropy_map == NULL) {
        fprintf(stderr, "Repeaty: Error: Failed to allocate memory for scrambled entropy map for sequence %s. Skipping.\n", header);
        free_interval_positions_map(original_interval_positions_map, MAX_INTERVAL_LEN, VALID_AA_COUNT);
        free_interval_counts(current_interval_counts, MAX_INTERVAL_LEN, VALID_AA_COUNT);
        free(sum_samelength_plog2p_by_length);
        free(total_intervals_at_length);
        return;
    }

    ScrambledEntropyStats scrambled_total_entropy_stats;
    ScrambledEntropyStats scrambled_same_entropy_stats;
    ScrambledEntropyStats scrambled_diff_entropy_stats;
    init_scrambled_entropy_stats(&scrambled_total_entropy_stats);
    init_scrambled_entropy_stats(&scrambled_same_entropy_stats);
    init_scrambled_entropy_stats(&scrambled_diff_entropy_stats);

    char scrambled_sequence_buf[MAX_SEQ_LEN + 1];
    Entropies current_scramble_entropies;

    for (int k = 0; k < num_scrambles_to_run; ++k) {
        if (verbose_flag && (k + 1) % 100 == 0) {
            fprintf(stderr, "    Processed %d/%d scrambles for %s.\n", k + 1, num_scrambles_to_run, header);
        }
        generate_scrambled_sequence(aa_counts, seq_len, scrambled_sequence_buf);

        // analyze_sequence_intervals still uses max_interval_len_for_seq_current for loop bounds
        int scramble_interval_count_val = analyze_sequence_intervals(scrambled_sequence_buf, seq_len, aa_counts, max_interval_len_for_seq_current, scrambled_entropy_map, &current_scramble_entropies, current_interval_counts);

        if (scramble_interval_count_val > 0) {
            if (!add_scrambled_plogp(&scrambled_total_entropy_stats, current_scramble_entropies.total_entropy) ||
                !add_scrambled_plogp(&scrambled_same_entropy_stats, current_scramble_entropies.same_residue_entropy) ||
                !add_scrambled_plogp(&scrambled_diff_entropy_stats, current_scramble_entropies.diff_residue_entropy)) {
                fprintf(stderr, "Repeaty: Error: Failed to add summary scrambled plogp value. Skipping remaining scrambles for %s.\n", header);
                break; // Exit scramble loop on error
            }
        } else { // If no intervals in this scramble, add 0 to maintain counts
            add_scrambled_plogp(&scrambled_total_entropy_stats, 0.0);
            add_scrambled_plogp(&scrambled_same_entropy_stats, 0.0);
            add_scrambled_plogp(&scrambled_diff_entropy_stats, 0.0);
        }
    }
    if (verbose_flag) {
        fprintf(stderr, "  Finished analyzing %d scrambled sequences for %s.\n", num_scrambles_to_run, header);
    }

    // --- 3. Prepare and Print Results ---
    IntervalData* all_intervals_for_seq_local = NULL;
    size_t num_intervals_for_seq_local = 0;
    size_t capacity_intervals_for_seq_local = 0;

    // Populate all_intervals_for_seq_local with data for sorting and printing
    if (original_interval_total_count > 0) {
        for (int len = 0; len <= max_interval_len_for_seq_current; len++) { // Loop up to current max
            for (int a_idx = 0; a_idx < VALID_AA_COUNT; a_idx++) {
                for (int y_idx = 0; y_idx < VALID_AA_COUNT; y_idx++) {
                    int original_count_for_type = original_interval_positions_map[len][a_idx][y_idx].count;

                    if (original_count_for_type > 0) {
                        double p = (double)original_count_for_type / original_interval_total_count;
                        double observed_plogp = -p * log2(p);

                        double p_same_length = 0.0;
                        double samelength_plog2p_val = 0.0;
                        if (total_intervals_at_length[len] > 0) {
                            p_same_length = (double)original_count_for_type / total_intervals_at_length[len];
                            if (p_same_length > 0) {
                                samelength_plog2p_val = -p_same_length * log2(p_same_length);
                            }
                        }

                        double p_value_for_graph_local = 1.0;
                        double mean_scrambled_plogp_local = 0.0;
                        double std_dev_scrambled_plogp_local = 0.0;
                        double z_score_local = 0.0;

                        ScrambledEntropyStats* ses = &scrambled_entropy_map[len][a_idx][y_idx];
                        if (ses->count > 0) { // Should be num_scrambles_to_run
                            mean_scrambled_plogp_local = ses->sum_plogp / ses->count;
                            
                            double variance = (ses->sum_sq_plogp / ses->count) - (mean_scrambled_plogp_local * mean_scrambled_plogp_local);
                            if (variance < 0 && fabs(variance) < 1e-9) variance = 0.0; 
                            else if (variance < 0) { 
                                fprintf(stderr, "Repeaty: Warning: Negative variance for mapped AA %d...%d length %d. Setting to 0.\n",
                                        a_idx, y_idx, len);
                                variance = 0.0;
                            }
                            std_dev_scrambled_plogp_local = sqrt(variance);

                            if (std_dev_scrambled_plogp_local < 1e-9) {
                                if (fabs(observed_plogp - mean_scrambled_plogp_local) < 1e-9) p_value_for_graph_local = 1.0;
                                else p_value_for_graph_local = 0.0;
                            } else {
                                z_score_local = (observed_plogp - mean_scrambled_plogp_local) / std_dev_scrambled_plogp_local;
                                double p_one_tail = normal_cdf(z_score_local);
                                double p_two_tail = (z_score_local < 0) ? (2 * p_one_tail) : (2 * (1.0 - p_one_tail));
                                if (p_two_tail > 1.0) p_two_tail = 1.0;
                                p_value_for_graph_local = (p_two_tail > 0.5) ? (1.0 - p_two_tail) : p_two_tail;
                            }
                        }

                        // Resize all_intervals_for_seq_local if needed
                        if (num_intervals_for_seq_local == capacity_intervals_for_seq_local) {
                            capacity_intervals_for_seq_local = (capacity_intervals_for_seq_local == 0) ? 10 : capacity_intervals_for_seq_local * 2;
                            IntervalData* temp = (IntervalData*)realloc(all_intervals_for_seq_local, capacity_intervals_for_seq_local * sizeof(IntervalData));
                            if (temp == NULL) {
                                fprintf(stderr, "Repeaty: Error: Failed to reallocate memory for interval data for sorting. Skipping remaining output for %s.\n", header);
                                free(all_intervals_for_seq_local); // Free partially allocated memory
                                all_intervals_for_seq_local = NULL;
                                num_intervals_for_seq_local = 0;
                                break; // Exit loop
                            }
                            all_intervals_for_seq_local = temp;
                        }

                        // Populate IntervalData struct
                        all_intervals_for_seq_local[num_intervals_for_seq_local].plogp_value = observed_plogp;
                        strncpy(all_intervals_for_seq_local[num_intervals_for_seq_local].seq_id, header, MAX_SEQUENCE_NAME_SIZE);
                        all_intervals_for_seq_local[num_intervals_for_seq_local].seq_id[MAX_SEQUENCE_NAME_SIZE] = '\0';
                        all_intervals_for_seq_local[num_intervals_for_seq_local].a_char = VALID_AMINO_ACIDS[a_idx];
                        all_intervals_for_seq_local[num_intervals_for_seq_local].y_char = VALID_AMINO_ACIDS[y_idx];
                        all_intervals_for_seq_local[num_intervals_for_seq_local].length = len;
                        all_intervals_for_seq_local[num_intervals_for_seq_local].count = original_count_for_type;
                        all_intervals_for_seq_local[num_intervals_for_seq_local].samelength_p = p_same_length;
                        all_intervals_for_seq_local[num_intervals_for_seq_local].samelength_plog2p = samelength_plog2p_val;

                        all_intervals_for_seq_local[num_intervals_for_seq_local].mean_scrambled_plogp = mean_scrambled_plogp_local;
                        all_intervals_for_seq_local[num_intervals_for_seq_local].std_dev_scrambled_plogp = std_dev_scrambled_plogp_local;
                        all_intervals_for_seq_local[num_intervals_for_seq_local].z_score = z_score_local;
                        all_intervals_for_seq_local[num_intervals_for_seq_local].p_value = p_value_for_graph_local;

                        num_intervals_for_seq_local++;
                    }
                }
            }
            if (all_intervals_for_seq_local == NULL && num_intervals_for_seq_local == 0) break; // Break outer loop if realloc failed
        }
    }

    // Declare and initialize significant_intervals_list at a wider scope
    IntervalData* significant_intervals_list = NULL;
    size_t num_significant_intervals_list = 0;
    size_t capacity_significant_intervals_list = 0;

    if (all_intervals_for_seq_local != NULL && num_intervals_for_seq_local > 0) {
        // Sort all intervals by P-value first to easily pick significant ones
        // sort_by_pvalue is now a fixed behavior, not an option
        qsort(all_intervals_for_seq_local, num_intervals_for_seq_local, sizeof(IntervalData), compareIntervalData);

        // Populate significant_intervals_list
        for (size_t i = 0; i < num_intervals_for_seq_local; i++) {
            if (all_intervals_for_seq_local[i].p_value <= SIGNIFICANCE_THRESHOLD &&
                all_intervals_for_seq_local[i].count >= min_count &&
                all_intervals_for_seq_local[i].length <= max_len_output) {

                // Add to significant_intervals_list
                if (num_significant_intervals_list == capacity_significant_intervals_list) {
                    capacity_significant_intervals_list = (capacity_significant_intervals_list == 0) ? 20 : capacity_significant_intervals_list * 2;
                    IntervalData* temp = (IntervalData*)realloc(significant_intervals_list, capacity_significant_intervals_list * sizeof(IntervalData));
                    if (temp == NULL) {
                        fprintf(stderr, "Repeaty: Error: Failed to reallocate memory for significant_intervals_list. Skipping.\n");
                        free(significant_intervals_list);
                        significant_intervals_list = NULL;
                        num_significant_intervals_list = 0;
                        break; // Exit loop
                    }
                    significant_intervals_list = temp;
                }
                significant_intervals_list[num_significant_intervals_list++] = all_intervals_for_seq_local[i];
            }
        }

        /* Short output format: print top 10 significant intervals by P-value and top 10 by count */
        /* --- Top 10 by P-value (already sorted from significant_intervals_list) --- */
        if (num_significant_intervals_list > 0) {
            fprintf(f4_out, "\nRepeaty: Top 10 Significant Intervals (by P-value) for %s:\n", header);
            print_short_interval_header(f4_out);

            for (size_t i = 0; i < num_significant_intervals_list && i < 10; i++) {
                 print_short_interval_line(f4_out, &significant_intervals_list[i], original_interval_positions_map);
            }
        } else {
            fprintf(f4_out, "\nRepeaty: No significant intervals to list by P-value for %s.\n", header);
        }

        /* --- Top 10 by Count (from significant intervals) --- */
        if (num_significant_intervals_list > 0) {
            qsort(significant_intervals_list, num_significant_intervals_list, sizeof(IntervalData), compareIntervalDataByCount);

            fprintf(f4_out, "\nRepeaty: Top 10 Significant Intervals (by Count) for %s:\n", header);
            print_short_interval_header(f4_out);

            for (size_t i = 0; i < num_significant_intervals_list && i < 10; i++) {
                print_short_interval_line(f4_out, &significant_intervals_list[i], original_interval_positions_map);
            }
        } else {
            fprintf(f4_out, "\nRepeaty: No significant intervals to list by Count for %s.\n", header);
        }
    } else {
        fprintf(f4_out, "Repeaty: No valid intervals found for sequence %s.\n", header);
    }

    // --- Print Summary Line ---
    double z_total = 0.0, p_total = 1.0;
    double z_same = 0.0, p_same = 1.0;
    double z_diff = 0.0, p_diff = 1.0;

    if (scrambled_total_entropy_stats.count > 0) {
        double mean_total = scrambled_total_entropy_stats.sum_plogp / scrambled_total_entropy_stats.count;
        double var_total = (scrambled_total_entropy_stats.sum_sq_plogp / scrambled_total_entropy_stats.count) - (mean_total * mean_total);
        if (var_total < 0 && fabs(var_total) < 1e-9) var_total = 0.0;
        double std_total = sqrt(var_total);
        if (std_total < 1e-9) {
            if (fabs(observed_total_entropy - mean_total) < 1e-9) z_total = 0.0;
            else z_total = (observed_total_entropy < mean_total) ? -INFINITY : INFINITY;
        } else {
            z_total = (observed_total_entropy - mean_total) / std_total;
        }
        double p_one_tail_total = normal_cdf(z_total);
        double p_two_tail_total = (z_total < 0) ? (2 * p_one_tail_total) : (2 * (1.0 - p_one_tail_total));
        if (p_two_tail_total > 1.0) p_two_tail_total = 1.0;
        p_total = (p_two_tail_total > 0.5) ? (1.0 - p_two_tail_total) : p_two_tail_total;
    }

    if (scrambled_same_entropy_stats.count > 0) {
        double mean_same = scrambled_same_entropy_stats.sum_plogp / scrambled_same_entropy_stats.count;
        double var_same = (scrambled_same_entropy_stats.sum_sq_plogp / scrambled_same_entropy_stats.count) - (mean_same * mean_same);
        if (var_same < 0 && fabs(var_same) < 1e-9) var_same = 0.0;
        double std_same = sqrt(var_same);
        if (std_same < 1e-9) {
            if (fabs(observed_same_residue_entropy - mean_same) < 1e-9) z_same = 0.0;
            else z_same = (observed_same_residue_entropy < mean_same) ? -INFINITY : INFINITY;
        } else {
            z_same = (observed_same_residue_entropy - mean_same) / std_same;
        }
        double p_one_tail_same = normal_cdf(z_same);
        double p_two_tail_same = (z_same < 0) ? (2 * p_one_tail_same) : (2 * (1.0 - p_one_tail_same));
        if (p_two_tail_same > 1.0) p_two_tail_same = 1.0;
        p_same = (p_two_tail_same > 0.5) ? (1.0 - p_one_tail_same) : p_one_tail_same; // Corrected two-tailed p-value for same
    }

    if (scrambled_diff_entropy_stats.count > 0) {
        double mean_diff = scrambled_diff_entropy_stats.sum_plogp / scrambled_diff_entropy_stats.count;
        double var_diff = (scrambled_diff_entropy_stats.sum_sq_plogp / scrambled_diff_entropy_stats.count) - (mean_diff * mean_diff);
        if (var_diff < 0 && fabs(var_diff) < 1e-9) var_diff = 0.0;
        double std_diff = sqrt(var_diff);
        if (std_diff < 1e-9) {
            if (fabs(observed_diff_residue_entropy - mean_diff) < 1e-9) z_diff = 0.0;
            else z_diff = (observed_diff_residue_entropy < mean_diff) ? -INFINITY : INFINITY;
        } else {
            z_diff = (observed_diff_residue_entropy - mean_diff) / std_diff;
        }
        double p_one_tail_diff = normal_cdf(z_diff);
        double p_two_tail_diff = (z_diff < 0) ? (2 * p_one_tail_diff) : (2 * (1.0 - p_one_tail_diff));
        if (p_two_tail_diff > 1.0) p_two_tail_diff = 1.0;
        p_diff = (p_two_tail_diff > 0.5) ? (1.0 - p_one_tail_diff) : p_one_tail_diff; // Corrected two-tailed p-value for diff
    }

    fprintf(f4_out, "\nRepeaty: Summary for %s: Total Entropy: %.2e (Z: %.2f, P: %.2e), Same Residue Entropy: %.2e (Z: %.2f, P: %.2e), Different Residue Entropy: %.2e (Z: %.2f, P: %.2e)\n",
           header, observed_total_entropy, z_total, p_total,
           observed_same_residue_entropy, z_same, p_same,
           observed_diff_residue_entropy, z_diff, p_diff);


    /* Free memory specific to this sequence's processing */
    // These now use MAX_INTERVAL_LEN
    free_interval_positions_map(original_interval_positions_map, MAX_INTERVAL_LEN, VALID_AA_COUNT);
    free_scrambled_entropy_map(scrambled_entropy_map, MAX_INTERVAL_LEN, VALID_AA_COUNT);
    
    // Free the summary scrambled entropy stats
    free(scrambled_total_entropy_stats.plogp_values);
    free(scrambled_same_entropy_stats.plogp_values);
    free(scrambled_diff_entropy_stats.plogp_values);

    free(all_intervals_for_seq_local);
    free(significant_intervals_list);
    free_interval_counts(current_interval_counts, MAX_INTERVAL_LEN, VALID_AA_COUNT);
    free(sum_samelength_plog2p_by_length);
    free(total_intervals_at_length);
}


int main(int argc, char* argv[]) {

    int c; /* For getopt */
    int errflg = 0, run_id=0; /* Error flag for getopt */
    int oop=0, opened=0;
    char prefix[MAX_FILE_NAME_SIZE+1];
    char outfile[MAX_FILE_NAME_SIZE+1];
    char tempopt[50]="", optlist[50]="";
    prefix[0] = '\0'; /* Initialize prefix to empty string */
    time_t t;
    char startrealtime[30];

    /* Initialize random seed */
    strcpy(startrealtime,__TIME__);
    srand(time(&t));

    /* Initialize VALID_AA_COUNT and the amino acid character to index mapping */
    VALID_AA_COUNT = strlen(VALID_AMINO_ACIDS);
    initialize_aa_mapping();

     while((c = getopt(argc, argv, "hvO:s:")) != -1) { // Removed 'm', 'x', 'L', 'l'
        switch(c) {
            case 'h':
                print_help(argv[0]);
                exit(0);
            case 'v':
                verbose = 1; sprintf(tempopt,"%c",optopt); strcat(optlist, tempopt);
                break;
            case 'O':
                oop=1;
                strncpy(prefix, optarg, MAX_FILE_NAME_SIZE);
                prefix[MAX_FILE_NAME_SIZE] = '\0'; /* Ensure null-termination */
                sprintf(tempopt,"%c",optopt); strcat(optlist, tempopt);
                break;
            case 's':
                num_scrambles_to_run = atoi(optarg);
                if (num_scrambles_to_run <= 0 || num_scrambles_to_run > N_SCRAMBLES_MAX) {
                    fprintf(stderr, "Repeaty: Error: Number of scrambles (-s) must be between 1 and %d.\n", N_SCRAMBLES_MAX);
                    errflg++;
                }
                sprintf(tempopt,"%c",optopt); strcat(optlist, tempopt);
                break;
            case ':': /* Option requires an argument but none was given */
                fprintf(stderr, "Repeaty: Error: option -%c requires a value\n", optopt); // Added "Repeaty: Error:"
                errflg++;
                break;
            case '?': /* Unrecognized option */
                fprintf(stderr, "Repeaty: Error: unrecognized option: -%c\n", optopt); // Added "Repeaty: Error:"
                errflg++;
                break;
            default: /* Should not happen */
                errflg++;
                break;
        }
    }

    /* If there were errors in option parsing, print help and exit */
    if (errflg) {
        print_help(argv[0]);
        return 1;
    }

    /* After getopt, optind points to the first non-option argument (the fasta file) */
    if (optind >= argc) {
        fprintf(stderr, "Repeaty: Error: No FASTA file provided.\n");
        print_help(argv[0]);
        return 1;
    }

    /* Declare and open the FASTA input file */
    FILE *fasta_file = fopen(argv[optind], "r");
    if (fasta_file == NULL) {
        fprintf(stderr, "Repeaty: Error: Could not open FASTA file '%s'.\n", argv[optind]);
        return 1;
    }


/* HANDLE OUTPUT FILE STREAM(S) */
if(strlen(optlist)==0) { strcpy(optlist,"default"); }
if(oop) { opened=0;
while(!opened) { run_id = (int) floor( rand()/100000); /* Generate a 5-digit random number */
                 snprintf(outfile, sizeof(outfile), "%s.%d.%s.repeaty.txt", prefix, run_id, optlist);
                 if(!file_exist(outfile)) { f4=fopen(outfile, "w"); opened=1; }
               } /* end of while !opened */
} /* end if oop */
else { f4 = stdout; }


    char line[MAX_LINE_LENGTH + 2];    /* Buffer to read lines, +2 for newline and null terminator */
    char header[MAX_SEQUENCE_NAME_SIZE + 1];  /* Buffer for sequence header */
    char sequence[MAX_SEQ_LEN + 1]; /* Buffer for protein sequence */
    int reading_sequence = 0;      /* Flag to indicate if we are currently reading a sequence */
    sequence[0] = '\0';            /* Initialize sequence string as empty */
    header[0] = '\0';              /* Initialize header string as empty */


    /* Loop through each line of the FASTA file */
    while (fgets(line, sizeof(line), fasta_file)) {
        /* Remove trailing newline character if present */
        line[strcspn(line, "\n")] = 0;

        char* trimmed_line = line;
        /* Skip any leading whitespace characters */
        while (isspace((unsigned char)*trimmed_line)) {
            trimmed_line++;
        }

        /* Check if the line is a FASTA header (starts with '>') */
        if (trimmed_line[0] == '>') {
            /* If we were previously reading a sequence, process it now */
            if (reading_sequence) {
                /* Only process if there was a sequence */
                if (strlen(sequence) > 0) {
                    process_sequence_data(header, sequence, f4, verbose, min_count_threshold, max_output_interval_len);
                }
            }
            /* Start reading a new sequence */
            reading_sequence = 1;
            strncpy(header, trimmed_line + 1, MAX_SEQUENCE_NAME_SIZE);
            header[MAX_SEQUENCE_NAME_SIZE] = '\0'; /* Ensure null-termination */
            char* space_pos = strchr(header, ' ');
            if (space_pos != NULL) {
                *space_pos = '\0'; /* Truncate header at first space */
            }
            sequence[0] = '\0'; /* Clear sequence for new one */
        } else if (reading_sequence) {
            /* Append line to current sequence */
            for (int k = 0; trimmed_line[k] != '\0'; k++) {
                if (is_valid_amino_acid(trimmed_line[k])) { /* Only add valid amino acids */
                    if (strlen(sequence) < MAX_SEQ_LEN) {
                        char temp[2];
                        temp[0] = toupper((unsigned char)trimmed_line[k]); /* Store as uppercase */
                        temp[1] = '\0';
                        strcat(sequence, temp);
                    } else {
                        fprintf(stderr, "Repeaty: Warning: Sequence exceeding MAX_SEQ_LEN. Truncating.\n");
                        break;
                    }
                }
            }
        }
    }

    /* After the loop, process the very last sequence if one was being read */
    if (reading_sequence) {
        if (strlen(sequence) > 0) { /* Only process if there was a sequence */
            process_sequence_data(header, sequence, f4, verbose, min_count_threshold, max_output_interval_len);
        }
    }

    fclose(fasta_file);
    if (oop) fclose(f4);

    return 0;
}


