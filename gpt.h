/*
 * gpt.h
 *
 *  Created on: 2012/11/02
 *      Author: yamasita
 */

/*   nndegdType : Type of calculation of distance to the pixel with the same edge direction            */
/*        0 Original                                                                                   */
/*        1 Use Scan Table for both images                                                             */
/*        2 Use table (nndegdTbl) for the nearest point of the same direction and Scan table (ScanTbl) */
/*        3 Use table (nndegdTbl) for the nearest point of the same direction (The other is skipped)   */
/*                                                                                             */
/*   gptType Type of GPT calculation                                        */
/*        0 Original                                                        */
/*        1 Use GPT calcuation templates (gH0, gH1x, gH1y)                  */
/*                                                                                             */
/*   isTrunc                                                                */
/*        0 Not truncate Gauss integral                                     */
/*        1 Truncate Gauss integral with length of its parameter and lTrunc */
/*                                                                                             */
/*   lTrunc Coefficient to truncate Gauss integral (refer to isTrunc)       */

/* Table for calculating GPT template */
typedef struct {
	int     leng;     /* Length of image (one dimensional) */
	int     intgLeng; /* Length for integral*/
	int    *st;       /* Starting point referred by Gauss integral */
	int    *gSt;      /* Starting point of Gauss function for Gauss integral */
	int    *gEnd;     /* End point of Gauss function for Gauss integral      */
	double  gamma;    /* Gamma (Parameter of Gauss function*/
	double *gauss;    /* Gauss function */
} GptTblLine;

/* Image Information */
typedef struct {
	int        *img;       /* Pixel values */
	int        *ang;       /* Edge directions */
	double     *can;       /* Canonical image */
	double     *weight;    /* Weight for image */
	double     *gH0;       /* H0  of Gpt calculation template */
	double     *gH1x;      /* H1x of Gpt calculation template */
	double     *gH1y;      /* H1y of Gpt calculation template */
	double     *nndegdTbl; /* Distance table of the same direction */
  int        *imgDiler;    /* Pixel values */
  int        *angDiler;    /* Edge directions */
  double     *canDiler;       /* Canonical image */
  double     *weightDiler;    /* Canonical image */
  double     *nndegdTblDiler; /* Distance table of the same direction */
} ImgInf;

/* Gpt calculation context */
typedef struct {
	/* Fixed parameters */
	int         nx;            /* # of horizontal pixels */
	int         ny;            /* # of vertical pixels   */
	int         cx;            /* center of horizontal pixels */
	int         cy;            /* center of vertical pixels */
	int         nGamma;        /* # of values of Gauss integral parameter */
	double     *gammaList;     /* List of values of Gauss integeral paprameter */
	int         nDir;          /* # of quantized directions  (= 8) */
	double      lTrunc;        /* Parameter to cut off Gauss integeral*/
	int         smooth;         /* Make input image smooth by differential morphology*/
	/* Dnn Context */
	NndegdContext  *nndegdContext;
	/* Values fixed by the above parameters */
	int         nxy;           /* # of total pixels */
	int         nxy2;          /* # of total pixels of extend image for 3x3 filter */
	int         nxyDir;        /* nxy x nDir */
	int         nDirNd;        /* # of directions with no direction (9 / 1 for match without direction) */
	int         nxyDirNd;      /* nxy x nDirNd */
	int         nxyDirNdGamma; /* nxy x nDirG x nGamma */
	int         nst;           /* # of scan table */
	int         dDiagSq;       /* Squared length of diagonal of image */
	double      dDiag;         /* length of diagonal of image */
	/* Tables independent of content of image */
	GptTblLine *xFull;         /* Table to calculate GPT table pf x direction without truncation */
	GptTblLine *yFull;         /* Table to calculate GPT table pf y direction without truncation */
	GptTblLine *xTrunc;        /* Table to calculate GPT table pf x direction with truncation */
	GptTblLine *yTrunc;        /* Table to calculate GPT table pf y direction with truncation */
	/* Tables for the original algorithm */
	double     *gDist;       /* Distance of original algorithm */
	double     *gwt;         /* Guass function */
	/* Work region */
	double     *gh0;         /* Table of h0 */
	double     *gh1;         /* Table of h1 */
	double     *gh2;         /* Table of h2 */
	double     gpt[9];       /* Current projection transform parameter T*/
	double     gptNew[9];    /* New projection transform parameter */
	double     iGpt[9];      /* Inverse projection transform */
	int        *atrImg;      /* Image */
	int        *atrAng;      /* Angle */
	double     *atrCan;      /* Canonical image */
	double     *atrWeight;   /* Canonical image */
	int        *atrImgNew;   /* New image */
	int        *invAtrImg;   /* inverse transform of transform Image */
	double     *invAtrCan;   /* inverse transform of transform Canonical image */
	/* Parameter for execution */
	int         gptType;      /* GAT/GPT*/
	int         gptMethod;    /* Type of GPT calculation */
	int         gptDir;       /* direction of GPT */
	int         evalType;     /* Recognition type (correlation or nndegd) */
	int         evalDnnType;  /* Dnn type for evaluation */
	int         isTrunc;      /* Flat to truncate Gauss integral */
	int         nDiler;       /* # of dilation/erosion list */
	int        *dilerParmL;   /* Parameter list of dilation/erosion */
	int         nDiffEr;       /* # of iterations for eroded difference evaluation */
	/* Data changed with execution */
	ImgInf     *inpImgInf;    /* Information of test image */
	ImgInf     *tgtImgInf;    /* Information of registered image */
	/* Inverse check */
	int         inverseChk;   /* 0: skip inverse check, 1 : dp inverse check*/
	/* Transform check */
	int         transChk;  /* 0: skip transform check, 1 : dp transform check */
	/* Output */
	int        iter;          /* # of iteration */
	double     nndegd;           /* The minimum distance of the same direction */
	double     smpEval;       /* Simple correlation */
	double     gptEval;       /* GPT correlation */
	int        smpDilerN;      /* The best dilation/erosion no. */
	int        gptDilerN;      /* The best dilation/erosion no. */
	/* Extension */
	int         nxExt;         /* # of horizontal pixels of extended image */
	int         nyExt;         /* # of vertical pixels of extended image  */
	int         nExt;          /* # of extension pixels (left, right, top, bottom) */
	/* Work Directory */
	char       *workDir;       /* Directory for work files */
	/* Work image */
	double *workImgD1, *workImgD2, *workImgD3, *workImgD4;
	/* For debug */
	int        *diffErImg;     /* Difference image (Error image) */
	int        debug;          /* Flag for debug */
} GptContext;

void gptCorInit(int nx, int ny, int nExt, int nGamma, double *gammaList, GptContext *gptContext);
int  calGptCor(GptContext *gptContext);
int  gptCor(int gptType, int *g_ang1, double *g_can1, int *g_ang2, double *g_can2, double nndegd, double *gpt, GptContext *gptContext);
void mkGptTbl(GptContext *gptContext);
int calIntgLeng(int leng, double gamma, double lTrunc);
void mkGptTblLine(GptTblLine *gptTblLine);
void mkGptTemplate(ImgInf *imgInf, GptContext *gptContext);

void gptTblPr(GptContext *context);
void gptTemplatePr(GptContext *context);

