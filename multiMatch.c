/*
 * multiMatch.c
 *
 *  Created on: 2012/11/02
 *      Author: Yukihiko Yamashita
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>

#include "parameters.h"
#include "eval.h"
#include "gpt.h"
#include "multiMatch.h"
#include "utility.h"

//#define DEBUG
//#define DEBUGKNN
//#define DEBUGSORT1
//#define DEBUGSORT2

/* Heap sort*/
void heapSort(SortRecord **recordP, int nData);
void mkHeap(SortRecord **recordP, int leaf, int root);

MultiMatchContext *multiMatchInit(char *execName) {
	int loop;
	int rg_digit;

	/* Allocate regions for mmContext and gptContext */
	MultiMatchContext *mmContext = (MultiMatchContext *) malloc(sizeof(MultiMatchContext));
	GptContext *gptContext       = (GptContext *) malloc(sizeof(GptContext));

	/* Parameter list of Gauss integral (Descending order) */
	double  gammaList[]   = GAMMALIST;

	/* # of registered and testimages for each category */
	int nRgInCat[NCAT]      = NRGPAT;
	int nTsInCat[NCAT]      = NTSPAT;

	/* Diration and erorsion */
	int dilerParm[]         = DILER;

	/* Initialize Multi-Match Context*/
	mmContext->execName   = execName;
	/* Set work directory */
	mmContext->workDir    = WORKDIR;
	/* Set GPT context */
	mmContext->gptContext = gptContext;
	/* Calculation of # of all patterns */
	mmContext->nRgTatal = 0; mmContext->nRgMaxCat = 0; mmContext->nTsTatal = 0;
	for (loop = 0 ; loop < NCAT ; ++loop) {
		mmContext->nRgInCat[loop] = nRgInCat[loop];
		mmContext->nTsInCat[loop] = nTsInCat[loop];
		mmContext->nRgTatal      += nRgInCat[loop];
		if (mmContext->nRgMaxCat < nRgInCat[loop]) mmContext->nRgMaxCat = nRgInCat[loop];
		mmContext->nTsTatal      += nTsInCat[loop];
	}
	mmContext->nRgMaxInBatch = NRGMAXINBATCH;
	mmContext->nBatchRgTatal = 0;
	for (rg_digit = 0 ; rg_digit < NCAT ; ++rg_digit) {
		mmContext->nBatchRgTatal += mmContext->nBatchRgInCat[rg_digit]
		                          = (mmContext->nRgInCat[rg_digit] - 1) / mmContext->nRgMaxInBatch + 1;
	}

	/* Initialize gptContext */
	/* Set GPT calculation type */
	gptContext->isTrunc   = GAUSSMETHOD;   /* Set truncation */
	gptContext->gptMethod = GPTMETHOD;     /* Set gptMethod */
	gptContext->gptDir    = GPTDIR;        /* Set gptDir */
	/* Set recognition type */
	gptContext->evalType    = EVALTYPE;       /* Set evalType (EVALCOR or EVALNNDEGD) */
	gptContext->evalDnnType = EVALNNDEGDTYPE; /* Set nndegd type for EVALNNDEGD) */
	/* Set GPT calculation type */
	gptContext->gptType   = GPTTYPE;       /* Set gptType (GPT or GPT)*/
	/* Set work directory */
	gptContext->workDir   = WORKDIR;
	/* Dilation/Erorsion */
	gptContext->nDiler     = NDILER;
	gptContext->dilerParmL = (int *) malloc(sizeof(int) * gptContext->nDiler);
	for (loop = 0 ; loop < gptContext->nDiler ; ++loop)	gptContext->dilerParmL[loop] = dilerParm[loop];
	/* Dnn Context */
	gptContext->nndegdContext = (NndegdContext *) malloc(sizeof(NndegdContext));
	/* Dnn type */
	gptContext->nndegdContext->nndegdType = NNDEGDTYPE;

	/* Allocate region for calculation of GPT */
	gptCorInit(COL + 2 * NEXT, ROW + 2 * NEXT, NEXT, NGAMMA, gammaList, gptContext);

	/* Make tables independent of content of image */
	/* Make ScanTbl */
	mkScanTbl(gptContext->nndegdContext);
	/* Make tables to calculate GPT templates*/
	mkGptTbl(gptContext);

	/* Allocate region for registered images */
	mmContext->rgImgInfL = (ImgInf *) malloc(sizeof(ImgInf) * mmContext->nRgMaxInBatch);
	for (loop = 0 ; loop < mmContext->nRgMaxInBatch ; ++loop) {
		allcImgInf((mmContext->rgImgInfL + loop), gptContext);
	}

	/* Information for test images */
	mmContext->tsImgInf = (ImgInf *) malloc(sizeof(ImgInf));
	allcImgInf(mmContext->tsImgInf, gptContext);
	return mmContext;
}
/* Make image information file for registered and test images*/
int mkImgInfFile(MultiMatchContext *mmContext) {
	GptContext *gptContext = mmContext->gptContext; /* GPT Context pointer             */
	char    tsFile[MAX_FILENAME];      /* File name of test images                     */
	char    rgFile[MAX_FILENAME];      /* File name of registered images               */
	char    imgInfFile[MAX_FILENAME];  /* File name of information of registered image */
	int     rgCatN, rgNinCat;
	int     tsCatN, tsNinCat;

	FILE    *imgInfFp;                  /* File pointer of information of image         */

	/* Make templates for registered images (Out of measurement of execution time) */
	time(&(mmContext->startTime));
	printf("Start making registered image information (imgInf) file \n");
	for (rgCatN = 0 ; rgCatN <= 9 ; ++rgCatN) {
		printf("Making registered image information file for digit = %d \n", rgCatN);
		/* Open file of information of registered images */
		sprintf(imgInfFile, "%s%s/%s_rg_imgInf_%d", mmContext->workDir,VERSION, mmContext->execName, rgCatN);
		imgInfFp = fopen(imgInfFile, "wb");
		for (rgNinCat = 1 ; rgNinCat <= mmContext->nRgInCat[rgCatN] ; ++rgNinCat) {
			/* Read registered image */
			sprintf(rgFile, "%s/%s/ldg%d_%04d_gray.pgm", mmContext->workDir, PGMDIR, rgCatN, rgNinCat);
			/* Use the table of the minimum distance of same direction and not truncate Gauss integral */
			loadMkImgInf(rgFile, mmContext->rgImgInfL, gptContext);

			/* Save information of registered image to file */
			saveImgInf(imgInfFp, mmContext->rgImgInfL, gptContext);
		}
		fclose(imgInfFp);
	}
	time(&(mmContext->currentTime));
	printf("Time for making registered image information file %d (sec)\n",
			(int) (mmContext->currentTime - mmContext->startTime));

	/* Start to measure time for test image transform */
	time(&(mmContext->startTime));

	/* Make templates for test data (in measurement of execution time) */
	printf("Start making test image information file (imgInf) \n");
	for (tsCatN = 0 ; tsCatN <= 9 ; ++tsCatN) {
		printf("Making test image information file for digit = %d \n", tsCatN);
		/* Open file of information of test images */
		sprintf(imgInfFile, "%s%s/%s_ts_imgInf_%d", mmContext->workDir, VERSION, mmContext->execName, tsCatN);
		imgInfFp = fopen(imgInfFile, "wb");
		for (tsNinCat = 1 ; tsNinCat <= mmContext->nTsInCat[tsCatN] ; ++tsNinCat) {
			/* Read test image */
			sprintf(tsFile, "%s/%s/tdg%d_%04d_gray.pgm", mmContext->workDir, PGMDIR, tsCatN, tsNinCat);
			/* Use the table of the minimum distance of same direction and not truncate Gauss integral */
			loadMkImgInf(tsFile, mmContext->tsImgInf, gptContext);

			/* Save information of registered image to file */
			saveImgInf(imgInfFp, mmContext->tsImgInf, gptContext);
		}
		fclose(imgInfFp);
	}
	time(&(mmContext->currentTime));
	printf("Time for making test image information %d (sec)\n",
			(int) (mmContext->currentTime - mmContext->startTime));
	return 0;
}

int multiMatch(int rg_digit_start, int rg_digit_end, int rg_batch_start, int rg_batch_end, MultiMatchContext *mmContext) {
	GptContext *gptContext = mmContext->gptContext; /* GPT Context pointer             */
	char    rgImgInfFile[MAX_FILENAME];  /* File name of information of registered image */
	char    tsImgInfFile[MAX_FILENAME];  /* File name of information of test image */
	char    resultFile[MAX_FILENAME];  /* File name of result data                     */
#ifdef DEBUG
	char    tmpFile[MAX_FILENAME];     /* File name for temporary data                 */
#endif
	FILE    *rgImgInfFp;                /* File pointer of information of registered image */
	FILE    *tsImgInfFp;                /* File pointer of information of test image       */
	FILE    *resultFp;                  /* File pointer of result data                     */

	/* Information in process */
	int     tsCatN;     /* Category of current test image */
	int     tsNinCat;   /* Serial # in the category of current test image */
	int     tsNtotal;   /* Serial # for all of current test image */
	int     rgCatN;     /* Category of current registered image */
	int     nRgBatch;   /* Number of pattern in each batch */
	int     rgBatchN;   /* Batch no. for registered image */
	int     rgNinCat;   /* Serial # in the category of current registered image */
	int     rgNinBatch; /* position in batch                       */
	int     loop;       /* Loop variable */
	int     gptStatus;  /* Status of gpt correlation */

	double  smpEval, gptEval;   /* Evaluation of GPT matching */
	int     smpDilerN, gptDilerN;

	SortRecord    *smpRecord;   /* Record for sort of simple correlation */
	SortRecord    *gptRecord;   /* Record for sort of GPT correlation */
	SortRecord    **smpRecordP; /* Pointer of Pointer of record for sort of simple correlation */
	SortRecord    **gptRecordP; /* Pointer of Pointer of record for sort of GPT correlation */
	TopRecord     *topRecord;  /* Record for top N correlations */

	/* Allocate region for sort */
	smpRecord  = (SortRecord *) malloc(sizeof(SortRecord) * mmContext->nRgMaxInBatch);
	gptRecord  = (SortRecord *) malloc(sizeof(SortRecord) * mmContext->nRgMaxInBatch);
	smpRecordP = (SortRecord **) malloc(sizeof(SortRecord *) * mmContext->nRgMaxInBatch);
	gptRecordP = (SortRecord **) malloc(sizeof(SortRecord *) * mmContext->nRgMaxInBatch);

	for (loop = 0 ; loop < mmContext->nRgMaxInBatch ; ++loop) {
		smpRecordP[loop] = & (smpRecord[loop]);
		gptRecordP[loop] = & (gptRecord[loop]);
	}

	/* Allocate region to keep top N data */
	topRecord = (TopRecord *) malloc(sizeof(TopRecord));

	/* Start time measurement */
	time(&(mmContext->startTime));

	/* Process for a test image */
	for (rgCatN = rg_digit_start ; rgCatN <= rg_digit_end ; ++rgCatN) {
		/* Open information of register image of category rg_digit */
		printf("Open registered image information file for registered image digit = %d \n", rgCatN);
		sprintf(rgImgInfFile, "%s%s/%s_rg_imgInf_%d", mmContext->workDir, VERSION, mmContext->execName, rgCatN);

		rgImgInfFp = fopen(rgImgInfFile, "rb");
		if (rgImgInfFp == 0) {
			printf("registered image information file cannot be opened : %s \n", rgImgInfFile);
			exit(1);
		}

		/* Read information file of registered images (for one category) and close */
		for (rgBatchN = rg_batch_start ; rgBatchN < (rg_batch_end >= 0 ? rg_batch_end : mmContext->nBatchRgInCat[rgCatN]) ; ++rgBatchN) {
			nRgBatch = (rgBatchN != mmContext->nBatchRgInCat[rgCatN] - 1) ? mmContext->nRgMaxInBatch
					: (mmContext->nRgInCat[rgCatN] - rgBatchN * mmContext->nRgMaxInBatch);

			printf("Load registered image information file for registered image digit = %d batch = %d\n", rgCatN, rgBatchN);
			for (rgNinBatch = 0 ; rgNinBatch < nRgBatch ; ++rgNinBatch)
				loadImgInf(rgImgInfFp, &(mmContext->rgImgInfL[rgNinBatch]), gptContext);

			/* Prepare result file */
			sprintf(resultFile,  "%s%s/%s_result_%d_%d", mmContext->workDir, VERSION, mmContext->execName, rgCatN, rgBatchN);
			/* Make result file */
			printf("Open result file: %s\n", resultFile);
			resultFp = fopen(resultFile, "wb");
			if (resultFp == 0) {
				printf("Result File cannot be opened 1 \n");
				exit(1);
			}
			fclose(resultFp);
			/* Open result file */
			resultFp = fopen(resultFile, "r+b");
			if (resultFp == 0) {
				printf("Result File cannot be opened 2 \n");
				exit(1);
			}

			/* The category and sample for test image are initialized  */
			tsNtotal       = 0;
			/* Loop for ts digit */
			for (tsCatN = 0 ; tsCatN <= 9 ; ++tsCatN) {

				/* Read information of test image of category ts_digit */
				sprintf(tsImgInfFile, "%s%s/%s_ts_imgInf_%d", mmContext->workDir, VERSION, mmContext->execName, tsCatN);
				tsImgInfFp = fopen(tsImgInfFile, "rb");
				if (tsImgInfFp == 0) {
					printf("Test image information File cannot be opened \n");
					exit(1);
				}
				for (tsNinCat = 1 ; tsNinCat <= mmContext->nTsInCat[tsCatN] ; ++tsNinCat) {

					/* Read information file of a test image and close it */
					loadImgInf(tsImgInfFp, mmContext->tsImgInf, gptContext);

					rgNinCat = rgBatchN * mmContext->nRgMaxInBatch;
					for (rgNinBatch = 0 ; rgNinBatch < nRgBatch ; ++rgNinBatch) {
						++rgNinCat;

						/* Set data for test image transformation */
						gptContext->inpImgInf = mmContext->tsImgInf;
						gptContext->tgtImgInf = &(mmContext->rgImgInfL[rgNinBatch]);

						/* Calculate GPT correlation by test image transformation */
						gptStatus = calGptCor(gptContext); /* nndegdType = 3, gptType = 1 */
						if (gptStatus != 0) {
							printf("GPT fail on test image transform tsCatN = %d, tsNinCat = %d, rgCatN = %d, rgNinCat = %d iter = %d\n",
									tsCatN, tsNinCat, rgCatN, rgNinCat, gptContext->iter);
							gptContext->gptEval = -1.0e+10;
						}

#ifdef DEBUG
						/* Display the result and write to result file for test image transformation */
						printf("TsTrans ts(%d,%d) rg(%d, %d; %d) : %3d iters smpEval = %f gptEval = %f\n",
								tsCatN, tsNinCat, rgCatN, rgNinCat, gptContext->gptDilerN, gptContext->iter, gptContext->smpEval, gptContext->gptEval);
						sprintf(tmpFile, "%soutImage/rgTransGPT.pgm", mmContext->workDir);
						save_image_file(tmpFile, gptContext->atrImg, mmContext->gptContext->nx, mmContext->gptContext->ny);
#endif
						/* For only test patterns transformation */
						smpEval   = gptContext->smpEval;
						gptEval   = gptContext->gptEval;
						smpDilerN = gptContext->smpDilerN;
						gptDilerN = gptContext->gptDilerN;

						if (gptContext->gptDir != GPTFORTEST) {
							/* Set data for registered image transformation */
							gptContext->inpImgInf = &(mmContext->rgImgInfL[rgNinBatch]);
							gptContext->tgtImgInf = mmContext->tsImgInf;

							/* Calculate GPT evaluation by registered image transformation */
							gptStatus = calGptCor(gptContext); /* nndegdType = 3, gptType = 1 */
							if (gptStatus != 0) {
								printf("GPT fail on registered image transform tsCatN = %d, tsNinCat = %d, rgCatN = %d, rgNinCat = %d iter = %d\n",
										tsCatN, tsNinCat, rgCatN, rgNinCat, gptContext->iter);
								gptContext->gptEval = -1.0e+10;
							}

#ifdef DEBUG
							/* Display the result and write to result file for registered image transformation */
							printf("RgTrans ts(%d,%d) rg(%d, %d; %d) : %3d iters smpEval = %f gptEval = %f\n",
									tsCatN, tsNinCat, rgCatN, rgNinCat, gptContext->gptDilerN, gptContext->iter, gptContext->smpEval, gptContext->gptEval);
							sprintf(tmpFile, "%soutImage/rgTransGPt.pgm", mmContext->workDir);
							save_image_file(tmpFile, gptContext->atrImg, mmContext->gptContext->nx, mmContext->gptContext->ny);
#endif

						/* Select the larger evaluation of test or registered images */
							if (smpEval < gptContext->smpEval) {
								smpEval   = gptContext->smpEval;
								smpDilerN = - gptContext->smpDilerN;
							}
#ifndef NOGPT
							if (gptEval < gptContext->gptEval) {
								gptEval   = gptContext->gptEval;
								gptDilerN = - gptContext->gptDilerN;
							}
#endif
						}

#ifdef DEBUG
						printf("RgTrans ts(%d,%d) rg(%d, %d) : smpEval = %f gptEval = %f\n", 	tsCatN, tsNinCat, rgCatN, rgNinCat, smpEval, gptEval);
#endif

						/* Save simple correlations */
						smpRecordP[rgNinBatch]->rgCatN   = rgCatN;
						smpRecordP[rgNinBatch]->rgDilerN = smpDilerN;
						smpRecordP[rgNinBatch]->rgNinCat = rgNinCat;
						smpRecordP[rgNinBatch]->eval     = smpEval;
						/* Save GPT evaluation */
#ifndef NOGPT
						gptRecordP[rgNinBatch]->rgCatN   = rgCatN;
						gptRecordP[rgNinBatch]->rgDilerN = gptDilerN;
						gptRecordP[rgNinBatch]->rgNinCat = rgNinCat;
						gptRecordP[rgNinBatch]->eval     = gptEval;
#endif
					}

#ifdef DEBUGSORT1
					/* Display the larger correlations of simple and GPT */
					for (loop = 0 ; loop < nRgBatch ; ++loop) {
						printf("rgCatN = %d, tsCatN = %1d, tsNinCat = %4d, Smp = %f (%4d, %4d, %4d), Gpt = %f  (%4d, %4d, %4d)\n",
								rgCatN, tsCatN, tsNinCat,
								smpRecordP[loop]->eval, smpRecordP[loop]->rgCatN, , smpRecordP[loop]->rgNinCat, smpRecordP[loop]->rgDilerN
								gptRecordP[loop]->eval, gptRecordP[loop]->rgCatN, , gptRecordP[loop]->rgNinCat, gptRecordP[loop]->rgDilerN);
					}
#endif
					/* Heap sort of results for a test image */
					heapSort(smpRecordP, nRgBatch);
#ifndef NOGPT
					heapSort(gptRecordP, nRgBatch);
#endif
					/* Extract the best N data and save */
					topRecord->tsCatN   = tsCatN;
					topRecord->tsNinCat = tsNinCat;
					topRecord->tsNtotal = tsNtotal;
					for (loop = 0 ; loop < MAXK ; ++loop) {
						topRecord->smpEval[loop]     = smpRecordP[loop]->eval;
						topRecord->smpRgCatN[loop]   = smpRecordP[loop]->rgCatN;
						topRecord->smpRgDilerN[loop] = smpRecordP[loop]->rgDilerN;
						topRecord->smpRgNinCat[loop] = smpRecordP[loop]->rgNinCat;
						topRecord->gptEval[loop]     = gptRecordP[loop]->eval;
						topRecord->gptRgCatN[loop]   = gptRecordP[loop]->rgCatN;
						topRecord->gptRgNinCat[loop] = gptRecordP[loop]->rgNinCat;
					}

#ifdef DEBUGSORT2
					/* Display the larger correlations of simple and GPT */
					for (loop = 0 ; loop < MAXK ; ++loop) {
						printf("rgCatN = %d, tsCatN = %1d, tsNinCat = %4d, maxSmp = %f (%4d, %4d, %4d), maxGpt = %f  (%4d, %4d, %4d)\n",
								rgCatN, tsCatN, tsNinCat,
								topRecord->smpEval[loop], topRecord->smpRgCatN[loop], topRecord->smpRgDilerN[loop], topRecord->smpRgNinCat[loop],
								topRecord->gptEval[loop], topRecord->gptRgCatN[loop], topRecord->gptRgDilerN[loop], topRecord->gptRgNinCat[loop]);
					}
#endif
					fseek(resultFp, (long) sizeof(TopRecord) * tsNtotal, SEEK_SET);
					fwrite(topRecord, sizeof(TopRecord), 1, resultFp);

					++tsNtotal;
				}

				/* Display total time up to now */
				time(&(mmContext->currentTime));
				printf("Calculation time = %d (sec) until tsCatN = %d rgCatN = %d rgBatchN = %d \n",
						(int) (mmContext->currentTime - mmContext->startTime), tsCatN, rgCatN, rgBatchN);
				/* Close information file of test image of category ts_digit */
				fclose(tsImgInfFp);
			}
			fclose(resultFp);
		}
		fclose(rgImgInfFp);
	}

	/* Display total time */
	time(&(mmContext->currentTime));
	mmContext->totalTime = (int) (mmContext->currentTime - mmContext->startTime);
	printf("Total time (rg_digit from %d to %d) = %d (sec)\n", rg_digit_start, rg_digit_end, (int) (mmContext->currentTime - mmContext->startTime));
	return 0;
}

int multiKNN(int command, MultiMatchContext *mmContext, int kChk) {
	FILE *fp;
	char fileName[256];
	int kNN;
	int tsCatN, rgCatN, tsrgCatN, rgBatN, rgBatchNinCat;
	int order, totalOrder;
	int smpHit[NCAT], gptHit[NCAT];
	int smpLabel[NCAT], gptLabel[NCAT];
	int nPat[NCAT + 1];
	int smp_correct[NCAT + 1], gpt_correct[NCAT + 1];
	int sampleN, totalSampleN, tsNinTotal, rgNinSort;
	double smpCorrectRate[NCAT + 1], gptCorrectRate[NCAT + 1];
	int     nTsTotal;                /* # of all test samples */
	//GptContext    *gptContext = mmContext->gptContext; /* Gpt context */
	double smpErrorRateK[MAXK], gptErrorRateK[MAXK];
	double smpSelfEval, smpOtherEval, gptSelfEval, gptOtherEval;

	double topSmpEval[NCAT * NCAT]; /* For averaged evalations */
	double topGptEval[NCAT * NCAT]; /* For averaged evaluations */
	int smpConfMat[NCAT * NCAT];    /* Confusion matrix for simple matching */
	int gptConfMat[NCAT * NCAT];    /* Confusion matrix for GPT matching */
	TopRecord     *topRecord;       /* Region for top N results of each category and each test image*/
	TopRgCBRecord *topRgCBRecord;   /* Region for top N results of all categories */
	TopRgCBRecord *topRgCBRecordL;  /* Region for all results */
	SortRecord    *sortRecordL;     /* Region for sort of all data. Used through (smp/gpt)SortRecordP  */
	SortRecord    **sortRecordPL;   /* Region for pointer of sort of all data. Used for (smp/gpt)SortRecordP  */
	SortRecord    **smpSortRecordP; /* Region for sort of simple correlation */
	SortRecord    **gptSortRecordP; /* Region for sort of GPT correlation */
	size_t        nReadItems;       /* # of items read from fiel */

	nTsTotal = mmContext->nTsTatal;

	/* Allocate regions */
	topRecord     = (TopRecord *) malloc(sizeof(TopRecord));
	topRgCBRecordL = (TopRgCBRecord *) malloc(sizeof(TopRgCBRecord) * nTsTotal);

	/* Allocate regions for sort of all data */
	sortRecordL   = (SortRecord *) malloc(sizeof(SortRecord) * nTsTotal * MAXK * mmContext->nBatchRgTatal * 2);
	sortRecordPL  = (SortRecord **) malloc(sizeof(SortRecord *) * nTsTotal * MAXK * mmContext->nBatchRgTatal * 2);
	totalSampleN = 0;
	for (tsNinTotal = 0 ; tsNinTotal < nTsTotal ; ++ tsNinTotal) {
		topRgCBRecordL[tsNinTotal].smpSortRecordP = &(sortRecordPL[ 2 * tsNinTotal      * MAXK * mmContext->nBatchRgTatal]);
		topRgCBRecordL[tsNinTotal].gptSortRecordP = &(sortRecordPL[(2 * tsNinTotal + 1) * MAXK * mmContext->nBatchRgTatal]);
		for (sampleN = 0 ; sampleN < MAXK * mmContext->nBatchRgTatal  ; ++sampleN) {
			topRgCBRecordL[tsNinTotal].smpSortRecordP[sampleN] = &(sortRecordL[totalSampleN++]);
		}
		for (sampleN = 0 ; sampleN < MAXK * mmContext->nBatchRgTatal  ; ++sampleN) {
			topRgCBRecordL[tsNinTotal].gptSortRecordP[sampleN] = &(sortRecordL[totalSampleN++]);
		}
	}

	/* Load result data */
	rgBatchNinCat = 0;
	for (rgCatN = 0 ; rgCatN < NCAT  ; ++rgCatN) { /* Category of registered image */
		for (rgBatN = 0 ; rgBatN < mmContext->nBatchRgInCat[rgCatN]  ; ++rgBatN) { /* Category of registered image */
			/* Open result file */
			sprintf(fileName,  "%s%s/%s_result_%d_%d", mmContext->workDir, VERSION, mmContext->execName, rgCatN, rgBatN);
			fp = fopen(fileName, "rb");
			if (fp == NULL) {
				printf("can't open %s\n", fileName);
				exit(1);
			}

			for (tsNinTotal = 0 ; tsNinTotal < nTsTotal ; ++ tsNinTotal) {
#ifdef DEBUGKNN
				printf("ts_sample = %d\n", tsNinTotal);
#endif
				/* read recognition results one by one */
				nReadItems = fread(topRecord, sizeof(TopRecord), 1, fp);
				if (nReadItems == 0) {
					printf("Number of input data is not enough. tsNinTotal = %d\n", tsNinTotal);
					exit(0);
				}

				/* write data */
#ifdef DEBUGKNN
				printf("%d %d %d | %d %d :", topRecord->tsCatN, topRecord->tsNtotal, topRecord->tsNinCat, rgCatN, rgBatN);
				for (order = 0 ; order < MAXK ; ++order) {
					printf("%6.3f %6.3f | ", topRecord->smpEval[order], topRecord->gptEval[order]);
				}
				printf("\n");
				if (topRecord->tsNtotal != tsNinTotal) {
					printf("Data is inconsistent for ts_sample = %d\n", tsNinTotal);
					exit(1);
				}
#endif

				if (rgBatchNinCat == 0) {
					topRgCBRecord = &(topRgCBRecordL[tsNinTotal]);
					topRgCBRecord->tsCatN   = topRecord->tsCatN;
					topRgCBRecord->tsNinCat = topRecord->tsNinCat;
					topRgCBRecord->tsNtotal = topRecord->tsNtotal;
				}

				/* Collect results for each test image*/
				smpSortRecordP = topRgCBRecordL[tsNinTotal].smpSortRecordP;
				gptSortRecordP = topRgCBRecordL[tsNinTotal].gptSortRecordP;
				for (order = 0 ; order < MAXK ; ++order) {
					totalOrder = MAXK * rgBatchNinCat + order;
					smpSortRecordP[totalOrder]->eval     = topRecord->smpEval[order];
					smpSortRecordP[totalOrder]->rgCatN   = topRecord->smpRgCatN[order];
					smpSortRecordP[totalOrder]->rgNinCat = topRecord->smpRgNinCat[order];
					smpSortRecordP[totalOrder]->rgDilerN = topRecord->smpRgDilerN[order];
					gptSortRecordP[totalOrder]->eval     = topRecord->gptEval[order];
					gptSortRecordP[totalOrder]->rgCatN   = topRecord->gptRgCatN[order];
					gptSortRecordP[totalOrder]->rgNinCat = topRecord->gptRgNinCat[order];
					gptSortRecordP[totalOrder]->rgDilerN = topRecord->gptRgDilerN[order];
				}
				// printf("rg_digit = %d ts_digit = %d   %f \n", rgCatN, topNRecord->ts_digit, smpSortRecordP[totalOrder]->eval);
			}
			++rgBatchNinCat;
			fclose(fp);
		}
	}

	/* Re-sort */
	for (tsNinTotal = 0 ; tsNinTotal < nTsTotal ; ++ tsNinTotal) {
#ifdef DEBUGSORT
		topRgCBRecord = &(topRgCBRecordL[tsNinTotal]);
		printf("ts:%d.%d : \n",	topRgCBRecord->tsCatN, topRgCBRecord->tsNinCat);
		smpSortRecordP = topRgCBRecordL[tsNinTotal].smpSortRecordP;
		gptSortRecordP = topRgCBRecordL[tsNinTotal].gptSortRecordP;
		rgBatchNinCat = 0;
		for (rgCatN = 0 ; rgCatN < NCAT  ; ++rgCatN) { /* Category of registered image */
			for (rgBatN = 0 ; rgBatN < mmContext->nBatchRgInCat[rgCatN]  ; ++rgBatN) {
				for (order = 0 ; order < MAXK * NCAT; ++order) {
					totalOrder = MAXK * rgBatchNinCat + order;
					printf("%d,%d:%f %d,%d:%f | ", smpSortRecordP[totalOrder]->rgCatN, smpSortRecordP[totalOrder]->rgNinCat, smpSortRecordP[totalOrder]->eval,
							gptSortRecordP[totalOrder]->rgCatN, gptSortRecordP[totalOrder]->rgNinCat, gptSortRecordP[totalOrder]->eval);
				}
				++rgBatchNinCat;
			}
		}
		printf("\n");
#endif

		heapSort(topRgCBRecordL[tsNinTotal].smpSortRecordP, MAXK * mmContext->nBatchRgTatal);
		heapSort(topRgCBRecordL[tsNinTotal].gptSortRecordP, MAXK * mmContext->nBatchRgTatal);

#ifdef DEBUGSORT
		topRgCBRecord = &(topRgCBRecordL[tsNinTotal]);
		smpSortRecordP = topRgCBRecordL[tsNinTotal].smpSortRecordP;
		gptSortRecordP = topRgCBRecordL[tsNinTotal].gptSortRecordP;
		printf("ts:%d.%d : ",	topRgCBRecord->tsCatN, topRgCBRecord->tsNinCat);
		for (order = 0 ; order < MAXK ; ++order) {
			printf("%d,%d:%f %d,%d:%f | ", smpSortRecordP[order]->rgCatN, smpSortRecordP[order]->rgNinCat, smpSortRecordP[order]->eval,
					gptSortRecordP[order]->rgCatN, gptSortRecordP[order]->rgNinCat, gptSortRecordP[order]->eval);
		}
		printf("\n");
#endif
	}

	switch (command) {
	case 10: /* Output kNN results */
		for (kNN = 1 ; kNN <= MAXK ; ++kNN) {
			/* initialization of arrays */
			for (tsCatN = 0; tsCatN < NCAT + 1 ; ++tsCatN) {
				nPat[tsCatN] = 0;
				smp_correct[tsCatN] = 0;
				gpt_correct[tsCatN] = 0;
			}
			/* Initialize confusion matrix */
			for (tsrgCatN = 0; tsrgCatN < NCAT * NCAT ; ++tsrgCatN) {
				smpConfMat[tsrgCatN] = 0;
				gptConfMat[tsrgCatN] = 0;
			}

			for (tsNinTotal = 0; tsNinTotal < nTsTotal ; ++tsNinTotal) {
				topRgCBRecord = &(topRgCBRecordL[tsNinTotal]);
				++nPat[topRgCBRecord->tsCatN];
				++nPat[NCAT]; /* Total number of test patterns */

				/* k-Nearest Neighbor */
				/* Initialize the hits counter */
				for (rgCatN = 0; rgCatN < NCAT ; ++rgCatN) {
					smpHit[rgCatN]   = 0; smpLabel[rgCatN] = rgCatN;
					gptHit[rgCatN]   = 0; gptLabel[rgCatN] = rgCatN;
				}
				/* Count the registered pattern for each category in k */
				smpSortRecordP = topRgCBRecordL[tsNinTotal].smpSortRecordP;
				gptSortRecordP = topRgCBRecordL[tsNinTotal].gptSortRecordP;
				for (order = 0; order < kNN; ++order) {
					++smpHit[smpSortRecordP[order]->rgCatN];
					++gptHit[gptSortRecordP[order]->rgCatN];
				}
				/* Sort by the number of each category in k */
				quicksort_int(smpHit, smpLabel, 0, 9);
				quicksort_int(gptHit, gptLabel, 0, 9);

				if (smpLabel[0] == topRgCBRecord->tsCatN) {
					++smp_correct[topRgCBRecord->tsCatN];
					++smp_correct[NCAT]; /* For total recognition rate */
				}
				if (gptLabel[0] == topRgCBRecord->tsCatN) {
					++gpt_correct[topRgCBRecord->tsCatN];
					++gpt_correct[NCAT]; /* For total recognition rate */
				}
				++smpConfMat[topRgCBRecord->tsCatN * NCAT + smpLabel[0]];
				++gptConfMat[topRgCBRecord->tsCatN * NCAT + gptLabel[0]];
			}

			for (tsCatN = 0; tsCatN < NCAT + 1 ; ++tsCatN) {
				smpCorrectRate[tsCatN] = -1.0; /* Set value if nPat = 0 */
				gptCorrectRate[tsCatN] = -1.0; /* Set value if nPat = 0 */
				if (nPat[tsCatN] > 0) {
					smpCorrectRate[tsCatN] = (double) smp_correct[tsCatN] / nPat[tsCatN] * 100.0;
					gptCorrectRate[tsCatN] = (double) gpt_correct[tsCatN] / nPat[tsCatN] * 100.0;
				}
			}
			smpErrorRateK[kNN - 1] = smpCorrectRate[NCAT];
			gptErrorRateK[kNN - 1] = gptCorrectRate[NCAT];

			printf("Recognition rate by kNN (k = %d)\n", kNN);

			printf("Simple %5.3f : ", smpCorrectRate[NCAT]);
			for (tsCatN = 0; tsCatN < NCAT ; ++tsCatN) {
				printf("%6.3f ", smpCorrectRate[tsCatN]);
			}
			printf("\nGPT    %5.3f : ", gptCorrectRate[NCAT]);
			for (tsCatN = 0; tsCatN < NCAT ; ++tsCatN) {
				printf("%6.3f ", gptCorrectRate[tsCatN]);
			}
			printf("\nConfusion matrix by kNN (k = %d)\n", kNN);
			printf("Simple : \n");
			printf(" t\\r     0     1     2     3     4     5     6     7     8     9    # of test pats");
			for (tsCatN = 0; tsCatN < NCAT ; ++tsCatN) {
				printf("\n %1d : ", tsCatN);
				for (rgCatN = 0; rgCatN < NCAT ; ++rgCatN) {
					printf("%5.3f ", (double) smpConfMat[tsCatN * NCAT + rgCatN] / nPat[tsCatN]);
				}
				printf("  %5d", nPat[tsCatN]);
			}
			printf("\nGPT    : \n");
			printf(" t\\r     0     1     2     3     4     5     6     7     8     9    # of test pats");
			for (tsCatN = 0; tsCatN < NCAT ; ++tsCatN) {
				printf("\n %1d : ", tsCatN);
				for (rgCatN = 0; rgCatN < NCAT ; ++rgCatN) {
					printf("%5.3f ", (double) gptConfMat[tsCatN * NCAT + rgCatN] / nPat[tsCatN]);
				}
				printf("  %5d", nPat[tsCatN]);
			}
			printf("\n");
		}
		printf("\n\n");
		for (kNN = 1 ; kNN <= MAXK ; ++kNN) printf("%4.2f\n", 100.0 - smpErrorRateK[kNN - 1]);
		printf("\n");
		for (kNN = 1 ; kNN <= MAXK ; ++kNN) printf("%4.2f\n", 100.0 - gptErrorRateK[kNN - 1]);
		break;
	case 11: 	/* Show top MAXK evaluations of results */
		for (tsNinTotal = 0 ; tsNinTotal < nTsTotal ; ++ tsNinTotal) {
			topRgCBRecord = &(topRgCBRecordL[tsNinTotal]);
			smpSortRecordP = topRgCBRecord->smpSortRecordP;
			gptSortRecordP = topRgCBRecord->gptSortRecordP;
			printf("Test digit = %2d, serial = %5d\n ", topRgCBRecord->tsCatN, topRgCBRecord->tsNinCat);
			if (topRgCBRecord->tsCatN != smpSortRecordP[0]->rgCatN) {
				printf("X");
			} else {
				printf(" ");
			}
			for (order = 0; order < MAXK; ++order) {
				printf("%1d ", smpSortRecordP[order]->rgCatN);
			}
			printf(" : ");
			for (order = 0; order < MAXK; ++order) {
				printf("%5.3f ", smpSortRecordP[order]->eval);
			}
			printf(" : ");
			for (order = 0; order < MAXK; ++order) {
				printf("%4d ",smpSortRecordP[order]->rgNinCat);
			}
			printf("\n ");
			if (topRgCBRecord->tsCatN != gptSortRecordP[0]->rgCatN) {
				printf("X");
			} else {
				printf(" ");
			}
			for (order = 0; order < MAXK; ++order) {
				printf("%1d ", gptSortRecordP[order]->rgCatN);
			}
			printf(" : ");
			for (order = 0; order < MAXK; ++order) {
				printf("%5.3f ", gptSortRecordP[order]->eval);
			}
			printf(" : ");
			for (order = 0; order < MAXK; ++order) {
				printf("%4d ",gptSortRecordP[order]->rgNinCat);
			}
			printf("\n");
		}
		break;
	case 12: /* Output averaged evaluations */
		/* Initialize for the calculation of averaged evaluations */
		for (tsCatN = 0 ; tsCatN < NCAT  ; ++ tsCatN) {
			nPat[tsCatN] = 0;
			for (rgCatN = 0 ; rgCatN < NCAT  ; ++ rgCatN) {
				topSmpEval[tsCatN + NCAT * rgCatN] = 0;
				topSmpEval[tsCatN + NCAT * rgCatN] = 0;
			}
		}

		/* Calculation of averaged correlation (CAT times CAT matrix) */
		for (tsNinTotal = 0 ; tsNinTotal < nTsTotal ; ++ tsNinTotal) {
			topRgCBRecord = &(topRgCBRecordL[tsNinTotal]);
			smpSortRecordP = topRgCBRecord->smpSortRecordP;
			gptSortRecordP = topRgCBRecord->gptSortRecordP;
			for (rgCatN = 0 ; rgCatN < NCAT  ; ++ rgCatN) {
				/* Search the most matched registered image of category rgCatN */
				for (rgNinSort = 0 ; rgNinSort <  MAXK * mmContext->nBatchRgTatal ; ++rgNinSort) {
					if (rgCatN == smpSortRecordP[rgNinSort]->rgCatN) break;
				}
				topSmpEval[topRgCBRecord->tsCatN + NCAT * rgCatN] += smpSortRecordP[rgNinSort]->eval;
				for (rgNinSort = 0 ; rgNinSort <  MAXK * mmContext->nBatchRgTatal ; ++rgNinSort) {
					if (rgCatN == gptSortRecordP[rgNinSort]->rgCatN) break;
				}
				topGptEval[topRgCBRecord->tsCatN + NCAT * rgCatN] += gptSortRecordP[rgNinSort]->eval;
			}
			++nPat[topRgCBRecord->tsCatN];
		}

		printf("Averaged top evaluation simple in each registered for each test \n");
		printf(" t\\r     0     1     2     3     4     5     6     7     8     9\n");
		for (tsCatN = 0 ; tsCatN < NCAT  ; ++ tsCatN) {
			printf(" %1d : ", tsCatN);
			for (rgCatN = 0 ; rgCatN < NCAT  ; ++ rgCatN) {
				printf("%4.3f ", topSmpEval[tsCatN + NCAT * rgCatN] / nPat[tsCatN]);
			}
			printf("\n");
		}
		printf("\nAveraged top evaluation GPT in each registered for each test \n");
		printf(" t\\r     0     1     2     3     4     5     6     7     8     9\n");
		for (tsCatN = 0 ; tsCatN < NCAT  ; ++ tsCatN) {
			printf(" %1d : ", tsCatN);
			for (rgCatN = 0 ; rgCatN < NCAT  ; ++ rgCatN) {
				printf("%4.3f ", topGptEval[tsCatN + NCAT * rgCatN] / nPat[tsCatN]);
			}
			printf("\n");
		}
		printf("\n");

		/* Calculation of averaged correlation Self and other */
		smpSelfEval = smpOtherEval = gptSelfEval = gptOtherEval = 0.0;
		for (tsNinTotal = 0 ; tsNinTotal < nTsTotal ; ++ tsNinTotal) {
			topRgCBRecord = &(topRgCBRecordL[tsNinTotal]);
			smpSortRecordP = topRgCBRecord->smpSortRecordP;
			gptSortRecordP = topRgCBRecord->gptSortRecordP;
			for (rgNinSort = 0 ; rgNinSort <  MAXK * mmContext->nBatchRgTatal ; ++rgNinSort) {
				if (topRgCBRecord->tsCatN == smpSortRecordP[rgNinSort]->rgCatN) break;
			}
			smpSelfEval += smpSortRecordP[rgNinSort]->eval;
			for (rgNinSort = 0 ; rgNinSort <  MAXK * mmContext->nBatchRgTatal ; ++rgNinSort) {
				if (topRgCBRecord->tsCatN != smpSortRecordP[rgNinSort]->rgCatN) break;
			}
			smpOtherEval += smpSortRecordP[rgNinSort]->eval;
			for (rgNinSort = 0 ; rgNinSort <  MAXK * mmContext->nBatchRgTatal ; ++rgNinSort) {
				if (topRgCBRecord->tsCatN == gptSortRecordP[rgNinSort]->rgCatN) break;
			}
			gptSelfEval += gptSortRecordP[rgNinSort]->eval;
			for (rgNinSort = 0 ; rgNinSort <  MAXK * mmContext->nBatchRgTatal ; ++rgNinSort) {
				if (topRgCBRecord->tsCatN != gptSortRecordP[rgNinSort]->rgCatN) break;
			}
			gptOtherEval += gptSortRecordP[rgNinSort]->eval;
		}
		printf("%6.3f %6.3f %6.3f %6.3f\n",
				smpSelfEval / nTsTotal, smpOtherEval / nTsTotal, gptSelfEval / nTsTotal, gptOtherEval / nTsTotal);
		break;
	case 14: /* Output serials if the top of all is not correct */
		for (tsNinTotal = 0 ; tsNinTotal < nTsTotal ; ++ tsNinTotal) {
			topRgCBRecord = &(topRgCBRecordL[tsNinTotal]);
			smpSortRecordP = topRgCBRecord->smpSortRecordP;
			gptSortRecordP = topRgCBRecord->gptSortRecordP;
			if (topRgCBRecord->tsCatN != gptSortRecordP[0]->rgCatN) { /* If answer is not correct */
				for (rgBatchNinCat = 1 ; rgBatchNinCat < mmContext->nRgMaxInBatch ; ++rgBatchNinCat) {
					if (topRgCBRecord->tsCatN == gptSortRecordP[rgBatchNinCat]->rgCatN) break;
				}
				printf("%1d %5d %1d %3d %5d %1d %3d %5d\n", topRgCBRecord->tsCatN, topRgCBRecord->tsNinCat,
						gptSortRecordP[0]->rgCatN, gptSortRecordP[0]->rgNinCat, gptSortRecordP[0]->rgDilerN,
						topRgCBRecord->tsCatN, gptSortRecordP[rgBatchNinCat]->rgNinCat, gptSortRecordP[rgBatchNinCat]->rgDilerN);
			}
		}
		break;
	case 15: /* Output serials if kNN is not correct (k = kChk) ./GPTv1 15 name k (./GPTv1 15 smnistDnnb 7) */
		for (tsNinTotal = 0; tsNinTotal < nTsTotal ; ++tsNinTotal) {
			topRgCBRecord = &(topRgCBRecordL[tsNinTotal]);
			/* k-Nearest Neighbor */
			/* Initialize the hits counter */
			for (rgCatN = 0; rgCatN < NCAT ; ++rgCatN) {
				gptHit[rgCatN]   = 0; gptLabel[rgCatN] = rgCatN;
			}
			/* Count the registered pattern for each category in k */
			gptSortRecordP = topRgCBRecordL[tsNinTotal].gptSortRecordP;
			for (order = 0 ; order < kChk; ++order) {
				++gptHit[gptSortRecordP[order]->rgCatN];
			}
			/* Sort by the number of each category in k */
			quicksort_int(gptHit, gptLabel, 0, 9);

			if (gptLabel[0] != topRgCBRecord->tsCatN) { /* if kNN is not incorrect */
				printf("%1d %5d %1d", topRgCBRecord->tsCatN, topRgCBRecord->tsNinCat, gptLabel[0]);
				for (rgNinSort = 0 ; rgNinSort < MAXK ; ++rgNinSort) {
					printf(" %1d %5d", gptSortRecordP[rgNinSort]->rgCatN, gptSortRecordP[rgNinSort]->rgNinCat);
				}
				printf("\n");
			}
		}
		break;
	}
	return 0;
}

/* Heap sort */
void heapSort(SortRecord **recordP, int nData) {
	int        root, leaf;
	SortRecord *tmpP;
	leaf = nData - 1;      /* Leaf (position is the index of array) */
	root = nData / 2 - 1;  /* Root*/

	/* Initial semi-ordered tree construction */
	while (root >= 0) {
		mkHeap(recordP, leaf, root);
		--root;
	}
	/* Extract the smallest element from semi-ordered tree */
	while (leaf > 0) {
		/* Exchange root and leaf */
		tmpP          = recordP[0];
		recordP[0]    = recordP[leaf];
		recordP[leaf] = tmpP;
		/* Recontruct semi-ordered tree */
		--leaf;
		mkHeap(recordP, leaf, 0);
	}
}

/* Make semi-ordered tree when root is not correct position */
void mkHeap(SortRecord **recordP, int leaf, int root) {
	SortRecord *tmpP;
	int min;
	int child = root * 2 + 1;
	while(child <= leaf) {
		if (child < leaf && recordP[child + 1]->eval < recordP[child]->eval) {
			min = child + 1;
		} else {
			min = child;
		}
		if (recordP[root]->eval <= recordP[min]->eval) break;
		tmpP          = recordP[root];
		recordP[root] = recordP[min];
		recordP[min]  = tmpP;
		root          = min;
		child         = root * 2 + 1;
	}
}
