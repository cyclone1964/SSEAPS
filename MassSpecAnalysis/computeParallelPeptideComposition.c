/*+C
 ******************************************************************
 * This program computes a peptide compositions given their Mass Spec
 * peak weights. These Masses are given in Daltons and are assumed
 * accurate to at least four decimal places.
 *
 * This is a simple brute force knapsack search that has been
 * parallelized usingposix threads to speed execution. It has been
 * designed to support web services by making the instantiations
 * uniquely identified by an ID, the first command line argument.
 * 
 * Usage: computePeptideComposition ID targetMass <targetMass>  ...
 *
 * where:
 *
 * ID Is a uniqe run ID. This is used to form the name of all internal
 * filenames according to the specification for the SSEAPS project.
 *
 * targetMass is the target mass of the peptide to 4 decimal places
 * <targetMass> are optional additional masses up to 8.
 *
 * This will print all the compositions. The output file is a .csv
 * file that can be read into R. It has a column for each acid type,
 * and that column holds the counts for that type for the found
 * composition. The last column holds the target mass
 *
 * OR
 *
 * computePeptideComposition
 *
 * with no command line arguments will cause the program to run
 * various internal test cases.
 *
 ******************************************************************
 */

/* Includes */
#include <math.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

/* For the threaded implementation */
#include <pthread.h>
#include <semaphore.h>

/* File-Scope Constants, Macros, and Enumerations */

/* This determines at what level of the recursion we STOP making threads */
#define THREAD_LEVEL (5)

/* 
 * This is the number of amino acid types that can be in a peptide:
 * note we have done away with IsoLeucine since it cannot be
 * distinguished from Leucine by it's mass 
 */
#define NUM_AMINO_ACID_TYPES (19)

/* 
 * And the maximum number of acids in a peptide: this caps the total
 * run time. 
 */
#define MAX_PEPTIDE_SIZE (20)

/* This is the tolerance in the mass to declare a match */
#define TOLERANCE (0.000)

/* File-Scope Type Definitions */

/* A data base entry for an amino acid, which is initialized below */
typedef struct {
  char *symbol;
  char *name;
  char *formula;
  double mass;
} AMINO_ACID_DATA;

/* 
 * This structure is used to pass arguments to the processType
 * recursive function. We use this so that we can wrap all the
 * arguments into one structure for creating threads using the
 * pthreads interface
 */
typedef struct {
  int numAcids;			/* The total number of acids so far */
  int typeIndex;		/* The index of the type */
  int typeCounts[NUM_AMINO_ACID_TYPES]; /* The current type count list */

  long currentMass;		/* The current Mass */
  long numCombinations;		/* Number of combinations attempted */
} TYPE_ARGUMENTS;

/* File-Scope Variables */

/* 
 * This is the number of matches: access semaphore protected along
 * with the file output 
 */
static int numMatches;

/* 
 * This is the maximum number of acids that are possible given the
 * input mass
 */
static int maxAcids;

/* 
 * This is the targetMass and tolerance converted to integers
 * read-only after once set. 
 */
static int tolerance;
static long targetMass;

/* This is used to keep the file prints from becoming intertwined */
static sem_t *printMutex;

/* This is the output file pointer, semaphore protected */
static FILE *outputFp = NULL;

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
static long typeMasses[NUM_AMINO_ACID_TYPES];

/* This is the output file type */

/* File-Scope Prototypes */
static void *processType(void *vTypeArgument);
static void printCounts(TYPE_ARGUMENTS *typeArgument);

/*+F
 ********************************************************
 * 
 * processType - try all possible number of contributions of a type 
 *
 * This function, which is invoked recursively, does all the work in
 * finding the peptide composition.. It assigns all possible values to
 * the input type that can be assigned and computes the current mass
 * from that and tests it. If it is too big, then it terminates the
 * loop and returns. If it is just right, it prints the counts and
 * returns. If it is too small, it invokes itself with the new current
 * mass and the counts properly reduced for the next type.
 *
 * In order to speed the process, this routine will spawn threads ot
 * execute lower levels of the recursion for the first few layers of
 * recursion. Otherwise, it recurses directly.
 *
 * In order to recurse in a thread-safe manner, this routine will
 * instantiate local variables to hold the inputs to lower levels of
 * recursion when threads are going to be generated. Otherwise it
 * operates on the input directly.
 *
 * Parameters:
 *
 * All the arguments are held in a single structure, with a pointer
 * passed in containing them. The specific members of that structure
 * are:
 *
 * numAcids - the number of amino acids already in the peptide
 * typeIndex - the type of Amino Acid processed at this level of recursion
 * currentMass - the currentMass, accumulated over prior types
 * numCombinations - the number of combinations tried so far (used for stats)
 * typeCounts - a vector of counts assigned for each type so far. 
 *
 * numAcids is thus the sum of typeCounts
 * 
 * Returns: NONE
 ********************************************************
 */
static void *processType(void *vTypeArguments)
{
  int typeIndex;
  int loopCount = 0;
  TYPE_ARGUMENTS *inputArguments = vTypeArguments;
  
  /* The base case: we have no types left to assign */
  if (inputArguments->typeIndex == NUM_AMINO_ACID_TYPES) return(NULL);

  /* 
   * Now, if the type index is lower than the specified threadLevel,
   * we implement this loop using threads for any recursions
   */
  if (inputArguments->typeIndex < THREAD_LEVEL) {

    /* Threat bookeeping information */
    int typeCount;
    int numThreads = 0;
    pthread_t threadIds[MAX_PEPTIDE_SIZE+1];

    /* To be thread-safe, we have to have arguments for each thread */
    TYPE_ARGUMENTS *threadArgument;
    TYPE_ARGUMENTS threadArguments[MAX_PEPTIDE_SIZE+1];

    /* Now try each type count possibility */
    for (typeCount=0, threadArgument = threadArguments;
	 typeCount<= maxAcids - inputArguments->numAcids;
	 typeCount++, threadArgument++) {

      /* 
       * Set the arguments for this call from the input 
       *
       * NOTE: Each time through the loop we are re-copying the input,
       * which we never modify.
       */
      *threadArgument = *inputArguments;
      
      /* 
       * Add this number of this acid to the structure. Note that the
       * post-increment on the index in the second call is NOT linked
       * to the looping variable
       */
      threadArgument->numAcids += typeCount;
      threadArgument->typeCounts[inputArguments->typeIndex] = typeCount;
      threadArgument->currentMass +=
	typeCount * typeMasses[threadArgument->typeIndex++];
      threadArgument->numCombinations = 1;
      
      /* If this mass is too big, then we are done with this loop */
      if (threadArgument->currentMass > targetMass+TOLERANCE) break;

      /* If we found a match, print it */
      if (threadArgument->currentMass >= targetMass-TOLERANCE) {
	printCounts(threadArgument);
	break;
      }

      /* Otherwise spawn a thread */
      pthread_create(&threadIds[numThreads++],
		     NULL,
		     processType,
		     threadArgument);
    }

    /* Join any threads we created and do bookeeping on the modified inputs */
    while(--numThreads >=0) {
      pthread_join(threadIds[numThreads],NULL);
      inputArguments->numCombinations +=
	threadArguments[numThreads].numCombinations;
    }
    
  } else {

    /* 
     * This bit of code is highly optimized to reduce operations to a
     * bare minimum given the limitations in the calling
     * convention. Specifically, declaring variables here slows it
     * down and this loop was re-factored a few times to make it
     * faster.
     */
    typeIndex = inputArguments->typeIndex++;
    while (inputArguments->numAcids <= maxAcids) {
      
      /* If this mass is too big, then we are done with this loop */
      if (inputArguments->currentMass > targetMass+TOLERANCE) break;
      
      /* If we found a match, print it */
      if (inputArguments->currentMass >= targetMass-TOLERANCE) {
	printCounts(inputArguments);
	break;
      }

      /* Otherwise let's process the next type Argument */
      processType(inputArguments);

      /* Add one of this type */
      loopCount++;
      inputArguments->numAcids++;
      inputArguments->numCombinations++;
      inputArguments->typeCounts[typeIndex]++;
      inputArguments->currentMass += typeMasses[typeIndex];
    }

    /* Put the input back the way we found it */
    inputArguments->typeIndex--;
    inputArguments->numAcids -= loopCount;
    inputArguments->currentMass -= loopCount * typeMasses[typeIndex];
    inputArguments->typeCounts[typeIndex] = 0;
  }

  return(NULL);
}
/*+F
 ********************************************************
 * 
 * Prints the counts (when so enabled) and increments the 
 * number of matches in a thread safe way.
 *
 * Parameters:
 *
 * TYPE_ARGUMENTS *typeArguments - the type arguments to be printed
 * 
 * Returns: NONE
 ********************************************************
 */
static void printCounts(TYPE_ARGUMENTS *typeArgument)
{
  int itype;

  /* Semaphore protect this */
  sem_wait(printMutex);

  numMatches++;
  if (outputFp != NULL) {
    for (itype=0;itype<NUM_AMINO_ACID_TYPES;itype++) 
      fprintf(outputFp,"%02d,",typeArgument->typeCounts[itype]);
    fprintf(outputFp,"%.4f\n",(double)(typeArgument->currentMass)/10000);
  }

  sem_post(printMutex);

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
  char *pName, *idName;
  char semName[128];
  char fileName[128];
  
  int itype,itry,index;
  
  float runTime;
  double inputMass;

  struct timeval startTime, endTime;

  TYPE_ARGUMENTS typeArguments;

  /* Load the initial arguments */
  typeArguments.numAcids = 0;
  typeArguments.typeIndex = 0;
  typeArguments.currentMass = 0;
  memset(typeArguments.typeCounts,0,sizeof(typeArguments.typeCounts));

  /* Set up the target mass table */
  for (index=0;index<NUM_AMINO_ACID_TYPES;index++)
    typeMasses[index] = round(aminoAcidData[index].mass * 10000);
  
  /* Parse the input arguments */
  if (argc > 1) {
    pName = argv[0]; argc--; argv++;
    if (argc < 1) {
      printf("USAGE: %s ID mass <mass> <mass> \n",pName);
      exit(1);
    }

    /* Create the semaphore */
    idName = argv[0]; argc--; argv++;
    sprintf(semName,"computePeptideCompositionMutex-%s",idName);
    printMutex = sem_open(semName,O_CREAT,777,1);

    itry = 0;
    while (argc > 0) {

      /* Get the new target mass */
      sscanf(argv[0],"%lf",&inputMass);
      targetMass = round(inputMass * 10000);

      /* From that compute the maximum number of acids */
      maxAcids = ceil(inputMass/aminoAcidData[0].mass);
      if (maxAcids > MAX_PEPTIDE_SIZE) {
	maxAcids = MAX_PEPTIDE_SIZE;
	printf("Clipping Peptide Length at %d\n",maxAcids);
      } else {
	printf("Max Peptide Length: %d\n",maxAcids);
      }

      /* Open the output file */
      sprintf(fileName,"Compositions-%s-%d.csv",idName,itry);
      if ((outputFp = fopen(fileName,"w")) == NULL) {
	printf("Unable to open file <%s>\n",fileName);
	exit(1);
      }

      /* Initialize as necessary */
      memset(&typeArguments,0,sizeof(typeArguments));

      /* Process the data */
      printf("Process weight %s\n",argv[0]);
      processType(&typeArguments);

      /* Close the file */
      fclose(outputFp);

      /* Next weight please */
      argc--; argv++; itry++;
    }
    /* Close the file */
    sem_close(printMutex);
    sem_unlink(semName);
    
    exit(0);
  }

  /* 
   * Instead, let's run for various lengths, making sure to have one and
   * only one match. By using the largest mass as the input insures 
   * the maximum number of tries.
   */
  printf("\n\n No command line arguments: run test cases .... \n\n");
  printf("Timing Numbers (Thread Level %d)\n\n",THREAD_LEVEL);
  printf(" #Acids RunTime\n");

  sem_unlink("testSem");
  printMutex = sem_open("testSem",O_CREAT,777,1);
  for (index = 3; index < 12; index++) {

    /* Set the maximum number of acids from the loop counter */
    maxAcids = MAX_PEPTIDE_SIZE;

    /* Initialize as necessary */
    sprintf(fileName,"TimingTestCase-%02d.csv",index);
    outputFp = fopen(fileName,"w");

    /* Initialize as necessary */
    memset(&typeArguments,0,sizeof(typeArguments));

    /* Set the target mass */
    targetMass = round(index * typeMasses[NUM_AMINO_ACID_TYPES-1]);
    
    gettimeofday(&startTime,NULL);
    processType(&typeArguments);
    gettimeofday(&endTime,NULL);

    fclose(outputFp);
    runTime = 1e-6*(endTime.tv_usec - startTime.tv_usec);
    runTime += endTime.tv_sec - startTime.tv_sec;
    printf(" %6d %7.3f (%d:%ld)\n",
	   index,runTime,numMatches,typeArguments.numCombinations);
  }
  printf("\n\n Test for redundancy\n\n");
  printf("     Mass #Matches #Combinations\n");

  for (itry=0;itry<8;itry++) {

    /* Let's set up a peptide with 14 acids randomly selected */
    inputMass = 0.0;
    maxAcids = 14;
    for (index=0;index<maxAcids; index++)
      inputMass += aminoAcidData[rand()%NUM_AMINO_ACID_TYPES].mass;

    /* Initialize as necessary */
    memset(&typeArguments,0,sizeof(typeArguments));

    sprintf(fileName,"RedundanceTestCase-%02d.csv",itry);
    outputFp = fopen(fileName,"w");

    /* The target mass */
    targetMass = round(inputMass * 10000);
    processType(&typeArguments);

    printf(" %.4f (%d:%ld)\n",
	   inputMass,numMatches,typeArguments.numCombinations);

    fclose(outputFp);
  }
  sem_close(printMutex);
  sem_unlink("testSem");
  printf("\n\nDone!\n");
}
