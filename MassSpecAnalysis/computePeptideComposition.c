/*+C
 ******************************************************************
 * This program computes a peptide compoistion give it's Mass Spec
 * Peak. This Mass Spec is given in ??'s and is assumed accurate to
 * four decimal places.
 *
 * This is a simple brute force knapsack search
 *
 * Usage: computePeptideComposition targetMass maxPeptides <OutputFile>
 *
 * where:
 *
 * targetMass is the target mass of the peptide to 4 decimal places
 * maxPeptides is the maximum number of peptides that can be combined. 
 * OutputFile is an optional file where it will print the compositions
 *
 * This will print all the compositions. The output file is a .csv
 * file that can be read into R. It has a column for each acid type,
 * and tha tcolumn holds the counts for that type for the found
 * composition.
 *
 * OR
 *
 * computePeptideComposition
 *
 * with no command line arguments will cause the program to run
 * various test cases.
 ******************************************************************
 */

/* Includes */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>

/* File-Scope Constants, Macros, and Enumerations */

/* This is the number of amino acid types that can be in a peptide */
#define NUM_AMINO_ACID_TYPES (19)

/* This is the tolerance in the mass to declare a match */
#define TOLERANCE (0.000)

/* File-Scope Type Definitions */
typedef struct {
  char *symbol;
  char *name;
  char *formula;
  double mass;
} AMINO_ACID_DATA;

/* File-Scope Variables */

/* This enables the printing of the sequences, which can be very verbose */
static int printFlag = 0;

/* These count statistics of the call */
static unsigned int numMatches = 0;
static unsigned long numCombinations = 0;

/* 
 * This is the maximum number of amino acids, which is the second of the
 * command line arguments
 */
static int maxAminoAcids=10;

/* This is the tolernace converted to integer */
static int tolerance;

static AMINO_ACID_DATA aminoAcidData[NUM_AMINO_ACID_TYPES] = {
  "G", "Glycine",        "C2H5NO2",       75.0669,
  "A", "Alanine",        "C3H7NO2",       89.0935,
  "S", "Serine",         "C3H7NO3",       105.0930,
  "P", "Proline",        "C5H9NO2",       115.1310,
  "V", "Valine",         "C5H11NO2",      117.1469,
  "T", "Threonine",      "C4H9NO3",       119.1197,
  "C", "Cysteine",       "C3H7NO2S",      121.1590,
  /*   "I", "Isoleucine",     "C6H13NO2",      131.1736, */
  "L", "Leucine",        "C6H13NO2",      131.1736,
  "N", "Asparagine",     "C4H8N2O3",      132.1184,
  "D", "Aspartate",      "C4H7NO4",       133.1032,
  "Q", "Glutamine",      "C5H10N2O3",     146.1451,
  "K", "Lysine",         "C6H14N2O2",     146.1882,
  "E", "Glutamate",      "C5H9NO4",       147.1299,
  "M", "Methionine",     "C5H11NO2S",     149.2124,
  "H", "Histidine",      "C6H9N3O2",      155.1552,
  "F", "Phenylalanine",  "C9H11NO2",      165.1900,
  "R", "Arginine",       "C6H14N4O2",     174.2017,
  "Y", "Tyrosine",       "C9H11NO3",      181.1894,
  "W", "Tryptophan",     "C11H12N2O2",    204.2262
};

/* This is integer versions of the weights */
static int typeMasses[NUM_AMINO_ACID_TYPES];

/* File-Scope Prototypes */
static void findPeptides(double inputMass, int maxAminoAcids,FILE *fp);
static void processType(int numLeft,
			int typeIndex,
			int targetMass,
			int currentMass,
			int *typeCounts,
			FILE *fp);

static void printCounts(int targetMass,int *typeCounts,FILE *fp);

/*+F
 ********************************************************
 * 
 * processType - try all possible number of contributions of a type 
 *
 * This function, which is invoked recursively, does all the work in
 * finding the peptides. It assigns all possible values to the input
 * type that can be assigned and computes the current mass from that and
 * tests it. If it is too big, then it terminates the loop and
 * returns. If it is just right, it prints the counts and returns. If
 * it is too small, it invokes itself with the new current mass and
 * the num left properly reduced for the next type.
 *
 * Parameters:
 *
 * int numLeft - number of peptides currently left to be assigned
 * int typeIndex - current type index
 * int currentMass - the current mass
 * int *typeCounts - the current count of the types (for printing)
 *
 * Returns: NONE
 ********************************************************
 */
static void processType(int numLeft,
			int typeIndex,
			int targetMass,
			int currentMass,
			int *typeCounts,
			FILE *fp)
{
  int typeCount;
  int newMass;

  /* The base case: we have no types left to assign */
  /* printf("Test Type (%d,%d,%d)\n",typeIndex,numLeft,currentMass); */
  if (typeIndex == NUM_AMINO_ACID_TYPES) return;

  for (typeCount=0;typeCount <= numLeft; typeCount++) {
    typeCounts[typeIndex] = typeCount;
    newMass = currentMass + typeCount * typeMasses[typeIndex];
    numCombinations++;
    if (newMass > targetMass+TOLERANCE) {
      typeCounts[typeIndex] = 0;
      return;
    }
    if (newMass >= targetMass-TOLERANCE) {
      printCounts(targetMass,typeCounts,fp);
      typeCounts[typeIndex] = 0;
      return;
    } else {
      processType(numLeft - typeCount,
		  typeIndex+1,
		  targetMass,
		  newMass,
		  typeCounts, fp);
    }
  }

  /* Reset the type counts */
  typeCounts[typeIndex] = 0;
}
/*+F
 ********************************************************
 * 
 * Prints the counts (when so enabled) and increments the 
 * number of matches.
 *
 * Parameters:
 *
 * int targetMass - the target Mass matched
 * int *typeCounts - tye count of each type of acid
 * FILE *fp - file to print composition into (NULL = none)
 *
 * Returns: NONE
 ********************************************************
 */
static void printCounts(int targetMass,int *typeCounts,FILE *fp)
{
  int itype;

  if (fp!= NULL) {
    for (itype=0;itype<NUM_AMINO_ACID_TYPES;itype++){
      fprintf(fp,"%02d",typeCounts[itype]);
      if (itype < NUM_AMINO_ACID_TYPES-1)
	fprintf(fp,",");
    }
    fprintf(fp,"\n");
  }

  /* Increment the number of matches */
  numMatches++;
}

/*+M
 ********************************************************
 *
 * This routine findes the amino acid composition for a peptide givin
 * it's mass in Daltons, limited to no more than maxAminoAcids
 * entries. 
 *
 * Parameters:
 *
 * double inputMass - input mass in Daltons
 * int maxAminoAcids - maximum number of acids in Peptide
 * FILE *fp - file to print compositions to (NULL = none)
 *
 * Returns: NONE
 *
 ********************************************************
 */
static void findPeptides(double inputMass,
			 int maxAminoAcids,
			 FILE *fp)
{
  int itype;
  int targetMass;
  int typeCounts[NUM_AMINO_ACID_TYPES];

  /* Convert input mass and tolerance to integers */
  targetMass = round(10000*inputMass);
  tolerance = round(TOLERANCE * 10000);

  /* 
   * Now we load a vector with the integer representation of the
   * masses and set the target mass as an integer.
   */
  for (itype=0;itype < NUM_AMINO_ACID_TYPES; itype++) 
    typeMasses[itype] = round(10000*aminoAcidData[itype].mass);

  /* Initialize the type counts to 0 */
  for (itype=0;itype<NUM_AMINO_ACID_TYPES;itype++)
    typeCounts[itype] = 0;

  /* 
   * Find the peptides, starting with 0 mass, 0 type Counts, and 0
   * mass so far.
   */
  processType(maxAminoAcids,0,targetMass,0,typeCounts,fp);
}
/* The main routine. It runs in two different modes:
 *
 * When invoked with NO command line arguments, it runs a bunch of 
 * test cases for timing purposes and prints out the run times.
 *
 * When invoked with command line arguments, there must be two and an
 * optional thirde
 *
 * Mass - mass to be matched, in daltons, interpreted to 4 decimal places
 * Limit - maximum number of acids in the peptide to be considered
 * <fileName> - name of an output file that will contain the compositions
 *
 */
int main(int argc, char**argv)
{
  int itype,itry;
  int maxAminoAcids;
  
  float runTime;

  double inputMass;

  FILE *fp = NULL;
  
  struct timeval startTime, endTime;
  
  /* Parse the input arguments */
  if (argc > 1) {
    if (argc < 3 ||
	sscanf(argv[1],"%lf",&inputMass) != 1 ||
	sscanf(argv[2],"%d",&maxAminoAcids) != 1) {
      printf("USAGE: %s mass maxCount\n",argv[0]);
      exit(1);
    }

    /* If an output file was specified. so be it */
    if (argc == 4) {
      if ((fp = fopen(argv[3],"w")) == NULL) {
	printf("Warning: Unable to open output file <%s>\n",argv[3]);
      } else {      
	printf(" Write Compositions to <%s>\n",argv[3]);
	for (itype=0;itype < NUM_AMINO_ACID_TYPES; itype++)
	  fprintf(fp,"%2s,",aminoAcidData[itype].symbol);
	fprintf(fp,"\n");
      }
    }
    
    /* Find the peptides for the input */
    printFlag = 1;
    findPeptides(inputMass,16,fp);
    printf(" Mass %.4lf has %u possible compositions out of %ld peptides\n",
	   inputMass,numMatches,numCombinations);
    exit(0);
  }

  /* 
   * Now, let's run for various lengths, making sure to have one and
   * only one match. By using the largest mass as the input insures 
   * the maximum number of tries.
   */
  printf("\n\n No command line arguments: run test cases .... \n\n");
  printf("Timing Numbers\n\n");
  printf(" #Acids RunTime\n");
  printFlag = 0;
  for (maxAminoAcids = 2; maxAminoAcids < 12 /* 16 */; maxAminoAcids++) {
    numMatches = numCombinations = 0;
    gettimeofday(&startTime,NULL);
    findPeptides(maxAminoAcids * aminoAcidData[NUM_AMINO_ACID_TYPES-1].mass,
		 20,NULL);
    gettimeofday(&endTime,NULL);
    runTime = 1e-6*(endTime.tv_usec - startTime.tv_usec);
    runTime += endTime.tv_sec - startTime.tv_sec;
    printf(" %6d %7.3f (%d:%ld)\n",maxAminoAcids,runTime,numMatches,numCombinations);
  }

  exit(1);
  printf("\n\n Test for redundancy\n\n");
  printf("     Mass #Matches #Combinations\n");

  for (itry=0;itry<8;itry++) {
    /* Let's set up a peptide with 14 acids randomly selected */
    inputMass = 0.0;
    for (maxAminoAcids=0;maxAminoAcids<14; maxAminoAcids++)
      inputMass += aminoAcidData[rand()%NUM_AMINO_ACID_TYPES].mass;

    numMatches = numCombinations = 0;
    findPeptides(inputMass,16,NULL);
    printf("%4.4f %8u %13ld\n",inputMass,numMatches,numCombinations);
  }
  printf("\n\nDone!\n");
}
