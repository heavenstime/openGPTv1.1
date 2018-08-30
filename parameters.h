/*
 * parameters.h
 *
 *  Created on: 2012/11/06
 *      Author: yamasita
 */

/* Work directory */
#define VERSION "v1.1"
#define WORKDIR "./gptWork/"

/*************************** Fixed parameter ********************************/
#define GPTORG    0      /* Original GAT/PPT algorithm                     */
#define GPTTMPLT  2      /* GAT/PPT algorithm with templates               */

#define GPTFORTESTREG 0  /* Test and registered images are transformed (larger is result) */
#define GPTFORTEST    1  /* Only test image is transformed */

#define NOGAUSSTRUNC 0  /* Gauss integral without truncation */
#define GAUSSTRUNC   1  /* Gauss integral with truncation    */

#define ALLGAT       0  /* Only GAT */
#define ALTGATPPT    1  /* GAT and PPT are done alternatively */
#define GATPPTAFTGAT 2  /* PPT after GAT */

/* Parameter mixing between GAT and PPT gptCor() */
#define GPTGAT     0  /* only GAT */
#define GPTPPT     1  /* only PPT */
#define GPTGATPPT  2  /* (GAT + PPT) / 2 */

/* Evaluation method */
#define EVALCOR              0 /* Recognition by correlation */
#define EVALNNDEGD           1 /* Recognition by nndegd */
#define EVALNNDEGDVAR        2 /* Recognition by variance of nndegd */
#define EVALNNDEGDSQR        3 /* Recognition by square root of nndegd */
#define EVALNNDEGDSQRT       4 /* Recognition by square root of nndegd */
#define EVALSSD              6 /* Sum of squared distance (L2 distance between patterns)*/
#define EVALDIFFER           7 /* Eroded difference of image */
#define EVALWDCH             8 /* Weighed direction code histgram */
#define EVALNNDEGDWEIGHT    11 /* NNDEGD with curvature weight */
#define EVALNNDEGDSQRWEIGHT 13 /* NNDEGD square with curvature weight  */

/* Input image/inf for matching in single match program*/
#define SMATCHIMAGE         -1 /* Input image */
#define SMATCHNUMERIMAGE     0 /* Input numerical image */
#define SMATCHNUMERINF       1 /* Input numerical information */
/********************************************************************************************/

/**********************************   Select method   **************************************/
//#define NOGPT                          /* Define only simple correlation without GPT */
#define DIRMATCH                        /* Define if match with direction */
#define GPTTYPE        GATPPTAFTGAT     /* ALLGAT, ATLGATPPT, GATPPTAFTGAT */
#define GPTDIR         GPTFORTESTREG    /* GPT direction : GPTFORTESTREG, GPTFORTEST */
#define GPTMETHOD      GPTTMPLT         /* GPTORG or GPTTMPLT. */
#define GAUSSMETHOD    NOGAUSSTRUNC     /* GAUSSTRUNC NOGUASSTRUNC */
#define NNDEGDTYPE     NNDEGDTMPLT      /* NNDEGDORG, NNDEGDTBL, NNDEGDTMPLT NNDEGDTMPLTHALF */
#define EVALTYPE       EVALNNDEGDWEIGHT /* EVALCOR, EVALNNDEGD EVALNNDEGDVAR EVALNNDEGDWEIGHT EVALNNDEGDSQR EVALDIFFER EVALWDCH */
#define EVALNNDEGDTYPE NNDEGDTMPLT      /* NNDEGDORG, NNDEGDTBL, NNDEGDTMPLT NNDEGDTMPLTHALF */
/********************************************************************************************/

/* Set iteration parameters */
#define TIMESFSTGAT   6     /* Times of first GAT for GATPPTAFTGAT */
#define MAX_ITER      35    /* # of iteration for GPT calc.            */
#define STOP_ITER     1.00  /* if corNew < STOP_ITER * col, stop loop   */

/* Parameter list of Gauss integral */
#define GAMMALIST {1/4.0, 1/8.0, 1/16.0, 1/32.0, 1/64.0}
#define NGAMMA    4

/* Number of iteration for eroded difference evaluation */
#define NDIFFER   0

/* Dilation/Erorsion */
#define NDILER 1                     /* The number of types of dilation/erosion (= 1 : no dilation/erosion) */
#define DILER  {0, -2, -4, -6, 2, 4} /* Since background is white, minus is dilation and plus is erosion , first parameter should be 0 */
#define DILCOEF 0.1                  /* Step size for differential equation of dilation/erosion */
#define DILCUTTHRESH -2.0            /* Match with dilated/eroded patterns when eval > DILCUTTHRESH */

/* Parameters for curvature weight of NNDEGD */
#define CURVSMOOTH    -2    /* # of iterations to smooth the image                */
#define CURVITER      7     /* # of iterations to smooth the covariance component */
#define CURVCOEFF     10.0

/* Inverse check to exculde the case when an image are truncated at boundary but correlation is high */
#define INVERSECHK    0     /* 0: skip inverse check, 1 : do inverse check */
#define SELFCORLIMIT  0.97  /* Check bound of image */

/* Chech transformation to excude the case the transformation is too large    */
#define TRANSCHK      0    /* 0: Skip transform check, 1 : do transform check */
#define TRANSLIMIT    0.5  /* Lower bound of det(AA^T) / tr(AA^T) */

/* Parameters for GPT */
#define NEXT          0       /* Extension pixels of images for input image extension  */
#define INPUTSMOOTH   0       /* Smooth parameter for all input data                   */
#define WGT           1.5     /* Coefficient for Gauss integral parameter              */
#define EPS           1.0e-8  /* value assumed to be zero                              */
#define DIRECTEPS     20.0    /* value assumed to be zero for robert8                  */
#define BLUERCOEFLIST {0.075114, 0.123841, 0.075114, 0.123841, 0.204180, 0.123841, 0.075114, 0.123841,  0.075114}
                               /* Coefficient list for 3x3 bluer for calculating weight */
#define WDCHGAUSS     {0.0,   0.009, 0.017, 0.009, 0.0, 0.009, 0.057, 0.105, 0.057, 0.009, 0.017, 0.105, 0.194, 0.105, 0.017, 0.009, 0.057, 0.105, 0.057, 0.009, 0.0,   0.009, 0.017, 0.009, 0.0}
#define WDCHDIRSMOOTH {0.2,   0.6,   0.2}

/* Select Data type (TESTIMAGE MNIST IPTP) */
#define MNIST /* TESTIMAGE IPTP MNIST */
/* Size of image single data   */
#ifdef TESTIMAGE
#define ROW 60          /* Vertical size of image    */
#define COL 40          /* Horizontal size of image  */
#define PGMDIR ""
#define NRGPAT {1}
#define NTSPAT {1}
#endif
#ifdef MNIST
#define ROW 28          /* Vertical size of image    */
#define COL 28          /* Horizontal size of image  */
#define PGMDIR "mnistPgm"
#define NRGPAT  {5923, 6742, 5958, 6131, 5842, 5421, 5918, 6265, 5851, 5949}
#define NTSPAT  { 980, 1135, 1032, 1010,  982,  892,  958, 1028,  974, 1009}
#endif
#ifdef MNISTCV
#define ROW 28          /* Vertical size of image    */
#define COL 28          /* Horizontal size of image  */
#define PGMDIR "mnistCvPgm"
#define NRGPAT  {4738, 5393, 4766, 4904, 4673, 4336, 4734, 5012, 4680, 4759}
#define NTSPAT  {1185, 1349, 1192, 1227, 1169, 1085, 1184, 1253, 1171, 1190}
#endif
#ifdef IPTP
#define ROW 24          /* Vertical size of image    */
#define COL 16          /* Horizontal size of image  */
#define PGMDIR "cdrom1g"
#define NRGPAT  {2535, 1691, 1658, 1736, 1188, 1996, 1977, 1410, 2366, 1428}
#define NTSPAT  {1485, 2392, 1939, 2043, 2403, 1641, 1831, 2056,  986, 1140}
#endif

/* Constant data */
#define NCAT           10          /* # of categories                           */
#define MAXK           10          /* # of data to be kept for kNN (max k)      */
#define BLACK          0           /* pixel value of black                      */
#define WHITE          255         /* pixel value of white                      */
#define DIRECTION      8           /* # of directions                           */
#define NoDIRECTION    8           /* Index for no direction pixel              */
#define NRGMAXINBATCH  1000        /* The maximum number of rg data for a batch */
#define MAX_IMAGESIZE  1024
#define MAX_BUFFERSIZE 256
#define MAX_FILENAME   256         /* Filename length limit                      */
#define PI             3.141592654 /* PI                                         */
#define GPTFAIL        10000000
#define NNDEGDFAIL     1.0e+10     /* Constant to show NNDEGD failure            */


/* Debug eroded difference */
//#define DEBUGDIFFER
//#define DEBUGWEIGHT

/* main for test */
//#define TESTMAIN
