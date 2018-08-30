/*
 * multiMatch.h
 *
 *  Created on: 2012/11/06
 *      Author: yamasita
 */


/* Top data for each input image */
typedef struct {
  int    tsCatN;               /* Category of test pattern */
  int    tsNinCat;             /* Pattern # in a category */
  int    tsNtotal;             /* Pattern # in whole patterns */
  double smpEval[MAXK];        /* Top MAXK evaluations without transformation */
  double gptEval[MAXK];        /* Top MAXK evaluations with transformation */
  int    smpRgCatN[MAXK];       /* Categories of top MAXK evaluations without transformation */
  int    smpRgDilerN[MAXK];     /* Dilation # of top MAXK evaluations without transformation */
  int    smpRgNinCat[MAXK];     /* Pattern # of top MAXK evaluations without transformation */
  int    gptRgCatN[MAXK];       /* Categories of top MAXK evaluations with transformation */
  int    gptRgDilerN[MAXK];     /* Dilation # of top MAXK evaluations with transformation */
  int    gptRgNinCat[MAXK];     /* Pattern # of top MAXK evaluations with transformation */
} TopRecord;

/* Date for sort */
typedef struct {
  double eval;     /* Evaluation */
  int    rgCatN;   /* Category of the registered image */
  int    rgDilerN; /* Dilation/erorsion no */
  int    rgNinCat; /* # in the category of the registered image */
} SortRecord;

/* Top N data for all test categories */
typedef struct {
	int    tsCatN;               /* Category of test pattern */
	int    tsNinCat;             /* Pattern # in a category */
	int    tsNtotal;             /* Pattern # in whole patterns */
	SortRecord **smpSortRecordP; /* Pointer for array to sort results without transformation for a test patter */
	SortRecord **gptSortRecordP; /* Pointer for array to sort results with transformation for a test patter */
} TopRgCBRecord;


typedef struct {
	GptContext *gptContext;      /* GPT context */
	int     nRgTatal;            /* # of all registered images                */
	int     nRgMaxCat;           /* Max # of registered images in categories  */
	int     nRgMaxInBatch;       /* Max # of registered images in a batch for process  */
	int     nBatchRgInCat[NCAT]; /* # of batches for each category */
	int     nBatchRgTatal;       /* # of total batches */
	int     nTsTatal;            /* # of all test images                      */
	int     nRgInCat[NCAT];      /* # of registered images for each category */
	int     nTsInCat[NCAT];      /* # of test images for each category */
	char *execName;      /* execution name */
	char *workDir;       /* Directory for work files */
	ImgInf *tsImgInf;    /* Information for test image */
	ImgInf *rgImgInfL;   /* Information for registered image */
	time_t  startTime;   /* Start time */
	time_t  currentTime; /* Current time */
	double  totalTime;  /* Total time */
} MultiMatchContext;

/* Initialize multi matching */
MultiMatchContext *multiMatchInit(char *execName);
/* Make files which store image information */
int mkImgInfFile(MultiMatchContext *mmContext);
/* Multi matching */
int multiMatch(int rg_digit_start, int rg_digit_end, int rg_batch_start, int rg_batch_end, MultiMatchContext *mmContext);
/* kNN */
int multiKNN(int command, MultiMatchContext *mmContext, int kChk);
