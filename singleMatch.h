/*
 * singleMail.h
 *
 *  Created on: 2012/11/15
 *      Author: yamasita
 */

typedef struct {
	int        type;              /* -1 : testImage,  0 : use image data, 1 : use image information */
	GptContext *gptContext;        /* GPT context */
	ImgInf     *rgImgInf;          /* Image information of registered image */
	ImgInf     *tsImgInf;          /* Image information of test image */
	int         rgCatN;            /* Category of registered image */
	int         rgNinCat;          /* Image number in a Category */
	int         tsCatN;            /* Category of test image */
	int         tsNinCat;          /* Image number in a Category*/
	int         nRgMaxInBatch;     /* Maxxmum number of image in a batch*/
	int        *dilerParmL;        /* Parameter list of dilation/erosion */
	int tsOrgImg[ROW * COL];       /* Original test image */
	int rgOrgImg[ROW * COL];       /* Original registered  image */
	int tsTransImg[ROW * COL];     /* Transformed image test image */
	int rgTransImg[ROW * COL];     /* Transformed test image */
	char tsFile[MAX_FILENAME];     /* File name of test images           */
	char rgFile[MAX_FILENAME];     /* File name of registered images     */
	char outImgFile[MAX_FILENAME]; /* File name of output image           */
	char *outFlag;                 /* Part of file name for 4 images */
	char execName[64];             /* Execution image */
#ifdef DEBUGDIFFER
	int  diffErImgTs[ROW * COL];   /* Eroded difference for test data transformation */
	int  diffErImgTr[ROW * COL];    /* Eroded difference for training data transformation */
#endif
	int weightInt[ROW * COL];       /* Evaluation weight in integer format */
} SingleContext;

/* Initialize single match */
SingleContext *singleMatchInit();
/* Single matching */
int singleMatch(SingleContext *sContext);
/* Save 4 images (Original register and test, transformed register and test)*/
int save4Image(SingleContext *sContext);
