/*
 * singleMail.c
 *
 *  Created on: 2012/11/02
 *      Author: yamasita
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include "eval.h"
#include "gpt.h"
#include "parameters.h"
#include "singleMatch.h"
#include "multiMatch.h"
#include "utility.h"

SingleContext *singleMatchInit() {
	SingleContext *sContext = (SingleContext *) malloc(sizeof(SingleContext));
	GptContext *gptContext  = (GptContext *) malloc(sizeof(GptContext));
	int loop;
	/* Diration and erorsion */
	int dilerParm[]     = DILER;

	/* Parameter list of Gauss integral (Descending order) */
	double  gammaList[] = GAMMALIST;

	/* Set GPT context */
	sContext->gptContext    = gptContext;
	sContext->nRgMaxInBatch = NRGMAXINBATCH;

	/* Initialize gptContext */
	/* Set work directory */
	gptContext->workDir = WORKDIR;
	/* Set GPT calculation type */
	gptContext->isTrunc   = GAUSSMETHOD;  /* Set truncation */
	gptContext->gptMethod = GPTMETHOD;   /* Set gptMethod */
	gptContext->gptDir    = GPTDIR;      /* Set gptDir */
	/* Set recognition type */
	gptContext->evalType    = EVALTYPE;       /* Set recogType (RECOGCOR or RECOGNNDEGD) (not used) */
	gptContext->evalDnnType = EVALNNDEGDTYPE; /* Set nndegd type for EVALNNDEGD) */
	/* Set GPT calculation type */
	gptContext->gptType = GPTTYPE;      /* Set gptType (GPT or GAT)*/
	/* Dnn Context */
	gptContext->nndegdContext = (NndegdContext *) malloc(sizeof(NndegdContext));
	/* Dnn type */
	gptContext->nndegdContext->nndegdType = NNDEGDTYPE;

	/* Dilation/Erorsion */
	gptContext->nDiler     = NDILER;
	gptContext->dilerParmL = (int *) malloc(sizeof(int) * gptContext->nDiler);
	for (loop = 0 ; loop < gptContext->nDiler ; ++loop)
		gptContext->dilerParmL[loop] = dilerParm[loop];

	/* Allocate region for calculation of GPT */
	gptCorInit(COL + 2 * NEXT, ROW + 2 * NEXT, NEXT, NGAMMA, gammaList, gptContext);
	//printf("single nxy = %d\n", gptContext->nxy);

	/* Make tables independent of content of image */
	/* Make ScanTbl */
	mkScanTbl(gptContext->nndegdContext);
	/* Make tables to calculate GPT templates*/
	mkGptTbl(gptContext);

	/* Information for test images (This should be after gptCorInit) */
	sContext->tsImgInf = (ImgInf *) malloc(sizeof(ImgInf));
	allcImgInf(sContext->tsImgInf, gptContext);

	/* Allocate region for registered images */
	sContext->rgImgInf = (ImgInf *) malloc(sizeof(ImgInf));
	allcImgInf(sContext->rgImgInf, gptContext);

	return sContext;
}

int singleMatch(SingleContext *sContext) {
	GptContext *gptContext = sContext->gptContext;
	ImgInf     *tsImgInf   = sContext->tsImgInf;
	ImgInf     *rgImgInf   = sContext->rgImgInf;
	int        nx          = gptContext->nx; /* # of horizontal pixels */
	int        ny          = gptContext->ny; /* # of vertual pixels    */
	int        nxy         = nx * ny;
	int        tsNinCat,  rgNinBatch;
	int     iterTs, iterRg, pos, gptStatus;
	double  gptEval, gptEvalTs, gptEvalRg;   /* Correlation of GPT matching */

	FILE    *tsImgInfFp,  *rgImgInfFp;
	time_t  startTime;   /* Start time */
	time_t  currentTime; /* Current time */
	/* Start to measure time for test image transform */
	time(&startTime);
	/* Load data */
	if (sContext->type <= 0) {
		/* Use the table of the minimum distance of same direction and not truncate Gauss integral */
		loadMkImgInf(sContext->tsFile, tsImgInf, gptContext);
		/* Use the table of the minimum distance of same direction and not truncate Gauss integral */
		loadMkImgInf(sContext->rgFile, rgImgInf, gptContext);
	} else {
		/* Read information file of a test image and close it */
		tsImgInfFp = fopen(sContext->tsFile, "rb");
		if (tsImgInfFp == 0) {
			printf("Test image information File cannot be opened : %s \n" , sContext->tsFile);
			exit(1);
		}
		for (tsNinCat = 1 ; tsNinCat <= sContext->tsNinCat ; ++tsNinCat) loadImgInf(tsImgInfFp, tsImgInf, gptContext);
		fclose(tsImgInfFp);

		rgImgInfFp = fopen(sContext->rgFile, "rb");
		if (rgImgInfFp == 0) {
			printf("registered image information file cannot be opened : %s \n", sContext->rgFile);
			exit(1);
		}
		for (rgNinBatch = 1 ; rgNinBatch <= sContext->rgNinCat ; ++rgNinBatch) loadImgInf(rgImgInfFp, rgImgInf, gptContext);
		fclose(rgImgInfFp);
	}

	/* Save original registered/test image */
	for (pos = 0 ; pos < nxy ; ++pos) {
		sContext->rgOrgImg[pos] = rgImgInf->img[pos];
		sContext->tsOrgImg[pos] = tsImgInf->img[pos];
	}

	//gptContext->debug = 1;
	/* Process for a test image */
	/**************************** Set data for test image transformation **************************/
	gptContext->inpImgInf = tsImgInf;
	gptContext->tgtImgInf = rgImgInf;
#ifdef DEBUGDIFFER
	gptContext->diffErImg = sContext->diffErImgTs; /* To keep eroded difference image */
#endif
	/* Calculate GPT correlation by test image transformation */
	gptStatus = calGptCor(gptContext);
	if (gptStatus != 0) {
		printf("GPT fail on test image transform \n");
	}
	gptEvalTs = gptContext->gptEval;
	iterTs    = gptContext->iter;

	/* Display result of test image transformation */
//	printf("Test image transform  %3d iters smpCor = %f gptCor = %f nndegd = %f \n",
//			gptContext->iter, gptContext->smpEval, gptContext->gptEval, gptContext->nndegd);

	/* Save transformed test image transformation */
	for (pos = 0 ; pos < nxy ; ++pos) sContext->tsTransImg[pos] = gptContext->atrImg[pos];

	/* Keep GPT correlation by test image transformation */
	gptEval = gptContext->gptEval;

	/**************************** Set data for registered image transformation **********************/
	gptContext->inpImgInf = rgImgInf;
	gptContext->tgtImgInf = tsImgInf;
#ifdef DEBUGDIFFER
	gptContext->diffErImg = sContext->diffErImgTr; /* To keep eroded difference image */
#endif
	/* Calculate GPT correlation by registered image transformation */
	gptStatus = calGptCor(gptContext);
	if (gptStatus != 0) {
		printf("GPT fail on registered image transform \n");
	}

	/* Display result of registered image transformation */
//	printf("Registered image transform  %3d iters smpCor = %f gptCor = %f nndegd = %f \n",
//			gptContext->iter, gptContext->smpEval, gptContext->gptEval, gptContext->nndegd);
	/* Save transformed registered image transformation */
	for (pos = 0 ; pos < nxy ; ++pos) sContext->rgTransImg[pos] = gptContext->atrImg[pos];

	/* Select the larger correlation of test or registered images */
	gptEvalRg = gptContext->gptEval;
	gptEval   = gptEvalRg > gptEvalTs ? gptEvalRg : gptEvalTs;
	iterRg    = gptContext->iter;

	/* Display the results */
	printf("smp = %f, gpt = %f  tsTrans gpt = %f (iter %d)  rgTrans gpt = %f (iter %d)\n",
			gptContext->smpEval, gptEval, gptEvalTs, iterTs, gptEvalRg, iterRg);

	/* Display total time */
	time(&currentTime);
	printf("Total time = %d (sec)\n", (int) (currentTime - startTime));
	return 0;
}

/* rgTran  tsOrg  rgOrg  tsTran*/
int save4Image(SingleContext *sContext) {
	GptContext *gptContext = sContext->gptContext;
	int        nx  = gptContext->nx;  /* # of horizontal pixels                    */
	int        ny  = gptContext->ny;  /* # of vertual pixels                       */
	int        img[(nx * 4 + 3) * ny];
	int        ix, iy, pos, posTsTrans, posRgOrg, posRgTrans, posTsOrg;
	for (pos = 0 ; pos < (nx * 4 + 3) * ny ; ++pos) img[pos] = 0;

	pos = 0;
	for (iy = 0 ; iy < ny ; ++iy) {
		posRgTrans = (nx * 4 + 3) * iy;
		posTsOrg   = posRgTrans + nx + 1;
		posRgOrg   = posTsOrg   + nx + 1;
		posTsTrans = posRgOrg   + nx + 1;
		for (ix = 0 ; ix < nx ; ++ix) {
			/*
			img[posRgOrg]   = sContext->rgImgInf->img[pos];
			img[posTsTrans] = sContext->rgImgInf->imgDiler[pos];
			img[posTsOrg]   = sContext->rgImgInf->imgDiler[pos + nx * ny];
			img[posRgTrans] = sContext->rgImgInf->imgDiler[pos + nx * ny * 2];
			*/
			img[posRgTrans] = sContext->rgTransImg[pos];
			img[posTsOrg]   = sContext->tsOrgImg[pos];
			img[posRgOrg]   = sContext->rgOrgImg[pos];
			img[posTsTrans] = sContext->tsTransImg[pos];
			++pos; ++posTsTrans; ++posRgOrg; ++posRgTrans; ++posTsOrg;
		}
	}
	//		save_image_file(fileName, sContext->tsTransImg, nx, ny);
	save_image_file(sContext->outImgFile, img, nx * 4 + 3, ny);
	return 0;
}
