
Patterny is a package of decipherment helpers to aid in hypothesis 
generation for intrinsically-disordered, low-complexity or 
compositionally-biased protein parts. These helpers are described 
in detail in the article: 

 Harrison, PM. "Patterny: A troupe of decipherment helpers for intrinsic 
 disorder, low complexity and compositional bias in proteins", Biomolecules, ... 


TO COMPILE: 
----------- 
Use the following command in the Patterny directory: 

  make 

...and this to remove compiled executables: 

  make clean 


TO RUN THE PATTERNY SCRIPT: 
--------------------------- 

Use: 

  ./patterny INPUT_FILE_IN_FASTA_FORMAT OUTPUT_PREFIX 

The script patterny is a BASH script. 
To make it a ZSH script, change the shebang at the top to !/bin/zsh . 
Similarly, for the sub-scripts bandy.sh and moduley.sh . 

The output appears inside the patterny directory. 


SCRIPT INPUT AND OPTIONS
------------------------ 

The INPUT_FILE_IN_FASTA_FORMAT is an input file in FASTA format (protein 
sequences). 

Please note that the maximum sequence size is 5000 residues. So, if 
you feed it Titin, it will only analyze the first 5000 residues. 

The OUTPUT_PREFIX is a prefix that prefixes each output file. 

The script has two options: 
  -help    prints help information 
  -cmodules <yes|no> ....whether to make and process Cmodules. Default is 'yes'


TEST EXAMPLE DATA
----------------- 

To test the script there are two small files of examples in the 
sub-directory Examples. Output files for these examples are also 
provided in Examples/output, to check against. There may be slight 
differences owing to different pseudo-random number generation. 

To run patterny on the test data, use: 

  ./patterny Examples/further-examples-disprot.fasta further-examples-disprot 

  ./patterny -cmodules no Examples/further-examples-cmodules-yeast.fasta further-examples-cmodules-yeast 

In the latter case, they are already CModules (i.e., compositional 
modules defined using Moduley), so the option '-cmodules no' is specified. 


READY-COMPILED EXECUTABLES
--------------------------  

In the subdirectory ./bin, there are subsubdirectories macosx and linux 
containing compiled versions for: Mac OS version 15.5, build version 24F74, 
and GNU/Linux 3.10.0-1160.80.1.el7.x86_64.  


OVERVIEW OF OUTPUT
------------------ 

There are two sorts of output: (i) for the input sequences, and (ii) for 
the CModules, 'compositional modules', or regions optimized for 
compositional bias. 

The generated output files are as follows: 

from moduley: 
 *.moduley.out.txt.            = lists of CModules and Boundary Sets 
 *.Cmodules.out.fasta          = CModule sequence tracts in FASTA format 

from bandy: 
 *.bandy.out.txt               = lists of compositional banding 

from blocky:
 *.blocky.out.txt              = blockiness calculations for input sequences
 *.Cmodules.blocky.out.txt     = blockiness calculations for CModules 

from repeaty: 
 *.repeaty.out.txt             = interval entropy (repetitiveness) calculations for input sequences
 *.Cmodules.repeaty.out.txt    = interval entropy (repetitiveness) calculations for CModules

from runny: 
 *.runny.out.txt               = homopeptide (hpep) calculations for input sequences 
 *.Cmodules.runny.out.txt      = homopeptide (hpep) calculations for CModules 


OUTPUT FORMATS 
--------------

MODULEY: 

# The output format *.moduley.out.txt is, for example: 

NAME/ACCESSION	(MODULEY)	CMODULES/BOUNDARY_SET	SEQLEN	START	END	PVALUE	SIGNATURE	SUBSEQUENCE 
:	:	:
disprot|DP00118r011	MODULEY	CMODULES	449	67	94	8.34e-05	{LQ}	LSILRHQNLLKELQDLALQGAKERTHQQ 
:	:	:
disprot|DP00118r011	MODULEY	BOUNDARY_SET	449	128	149	1.55e-04	{D}	DAAEKRDDFKEVEKSDEDSDGD 

The columns are: 
 $1 = sequence name/accession 
 $2 = program name
 $3 = type of annotation (either CMODULE or BOUNDARY_SET) 
 $4 = sequence length (SEQLEN)
 $5 = start of subsequence (START) 
 $6 = end of subsequence (END) 
 $7 = binomial P-value for the bias formed by the residues making up the ... 
 $8 = ... bias signature {in curly brackets} (SIGNATURE) 
 $9 = subsequence of the CMODULE or BOUNDARY_SET (SUBSEQUENCE)


# The output format *.Cmodules.out.fasta is, for example: 

:	:	:
>disprot|DP00118r011_{R}_362_391_1.24e-05
RSMRLSFRARGYGFRGPGLQLRRGWRPNSR

>disprot|DP00118r011_{LQ}_67_94_8.34e-05
LSILRHQNLLKELQDLALQGAKERTHQQ
:	:	:

The ">" line contains a compound name made from: 
>NAME/ACCESSION_SIGNATURE_START_END_PVALUE 


--------------------------------------------------------------------------------------------------

BANDY: 

# The output format *.bandy.out.txt is, for example: 

OUTPUT_CRITERION	(BANDY)	NAME/ACCESSION	{PRIMARY_BIAS}	SEQLEN	TRACT_NUMBER	PRIMARY_BIAS_COUNT	PRIMARY_BIAS_PVALUE	DPB	ZSCORE	NORMALP	BAND_COUNT	OUTLIER:_INTERVAL_LIST	BAND:_LIST
:	:	:
HIGHEST_BAND_NUMBER:	BANDY	Q9UBG3	{Q}	398	1	66	2.39e-24	194	-0.46	6.43e-01	6	OUTLIERS:	|	BANDS:	|50-68|81-104|144-169|184-229|271-300|280-296|
HIGHEST_ZSCORE:	BANDY	Q9UBG3	{Q}	398	1	66	2.39e-24	162	0.04	9.70e-01	3	OUTLIERS:	N/A	BANDS:	|50-68|85-104|285-300|
LOWEST_ZSCORE:	BANDY	Q9UBG3	{Q}	398	1	66	2.39e-24	4	-2.37	1.77e-02	2	OUTLIERS:	N/A	BANDS:	|85-104|285-300|
:	:	:

The columns are: 
 $1 = output criterion (there are three HIGHEST_BAND_NUMBER, HIGHEST_ZSCORE, LOWEST_ZSCORE)  
 $2 = program name 
 $3 = sequence name/accession 
 $4 = {PRIMARY_BIAS} 
 $5 = sequence length (SEQLEN)  
 $6 = TRACT_NUMBER (this is per input sequence) 
 $7 = PRIMARY_BIAS_COUNT, the number of residues for the primary-bias amino acid 
 $8 = PRIMARY_BIAS_PVALUE, binomial P-value for the primary-bias amino acid across the input sequence 
 $9 = DPB, the distance to perfect banding 
$10 = ZSCORE of the DPB relative to population of 1000 scrambled sets of bands 
$11 = NORMALP, from the ZSCORE
$12 = BAND_COUNT, number of bands 
$13+: 
   OUTLIER:_INTERVAL LIST, outlier intervals between bands removed using median absolute deviation
   BAND:_LIST, list of bands (delimited by '|', start-end); 


--------------------------------------------------------------------------------------------------

BLOCKY: 

# The output format *.blocky.out.txt is, for example: 

NAME/ACCESSION	SEQLEN	B(BLOCKINESS)	NORMALP,ZSCORE,(MEAN_B +/- STD)	BLOCKINESS_FOR_AATYPES:|AA,COUNT,BINOMIALP,B(BLOCKINESS)|... 
:	:	:
disprot|DP00065r001	555	0.07	6.39e-02,-1.85,(0.08+/-0.00)	|T,76,2.4e-24,0.239|G,64,7.4e-172,0.148|D,61,1.8e-35,0.215|K,55,3.7e-03,0.161|... 
:	:	:

The columns are: 
 $1 = sequence name/accession 
 $2 = sequence length (SEQLEN)
 $3 = blockiness value (B)  
 $4 = NORMALP,ZSCORE,(MEAN_B +/- STD) : these are relative to population of 1000 scrambled sequences 
 $5 = BLOCKINESS_FOR_AATYPES, a list delimited by '|' in the form: |residue_type,count_of_residues,binomial_P_value_for_bias,residue-specific_blockiness|... 


--------------------------------------------------------------------------------------------------

REPEATY: 

# The output format *.repeaty.out.txt is, for example: 

:	:	:
Repeaty: Top 10 Significant Intervals (by P-value) for sp|P36157|MSA2_YEAST:
SeqID	Interval	Positions	Length	Count	-plog2p	ScrambleMean	ScrambleStdDev	Z-score	P-value
sp|P36157|MSA2_YEAST	P...E	99-196,201-298,234-331,251-348,266-363	96	5	1.35e-03	3.82e-04	1.30e-04	7.42	1.19e-13
sp|P36157|MSA2_YEAST	H...Q	277-284,295-302,325-332	6	3	8.55e-04	3.35e-04	7.48e-05	6.95	3.76e-12
:	:	:
Repeaty: Top 10 Significant Intervals (by Count) for sp|P36157|MSA2_YEAST:
SeqID	Interval	Positions	Length	Count	-plog2p	ScrambleMean	ScrambleStdDev	Z-score	P-value
:	:	:
Repeaty: Summary for sp|P36157|MSA2_YEAST: Total Entropy: 1.45e+01 (Z: -3.06, P: 2.19e-03), Same Residue Entropy: 1.06e+00 (Z: 1.79, P: 7.29e-02), Different Residue Entropy: 1.34e+01 (Z: -2.56, P: 1.06e-02)


There are two lists of top 10 significant intervals either sorted on P-value 
or on the total count of intervals of that type. 
Each line in these lists details an interval that occurs >=3 times, its positions 
and its statistics relative to 1000 scrambled sequences. 
'-plog2p' is the interval entropy term for the specific interval. 

The Summary lines contain information about the overall Interval Entropy. 
There are values for 'Total Entropy' amd its significance relative 
to the scrambled sequences, 'Same Residue Entropy', which is 
the same calculation but just considering intervals between residues of the same 
type, and also similaryly a 'Different Residue Entropy'. 

Significant P-values with negative Z-scores indicates REPETITIVENESS. 



