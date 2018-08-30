/*
 * gptMain.c
 *
 *  Created on: 2012/11/17
 *      Author: yamasita
 */

/* gptMain.c                                                                                                   */
/* kNN by GPT correlation                                                                                      */
/* Obtain the Affine transform to maximize canonical correlation on ROW x COL image                            */
/*                                                                                                             */
/* How to use                                                                                                  */
/*  ./gptMain comamnd execName errFileName                                                                     */
/*     command = 21 : Match by a registration images and a test images described in list in file               */
/*  ./gptMain comamnd                                                                                          */
/*     command = 20 : Match by a registration image and a test image                                           */
/*  ./gptMain comamnd execNme                                                                                  */
/*     comamnd =  0 : Match using all data  (Use after command 5)                                              */
/*     command =  5 : Make image information file for registered and test images                               */
/*     command = 10 : Recognition by k-nearest neighbor                                                        */
/*     command = 11 : Show best 10 correlations of registered pattern for each test pattern                    */
/*     command = 12 : Show average of best 10 correlations of results                                          */
/*     command = 14 : Make list of test pattern of which nearest pattern is in incorrect category              */
/*  ./gptMain comamnd  execName start_rg_digit end_rg_digit                                                    */
/*     comamnd =  1 : Match rg_digits from start_rg_digit to end_rg_digit  (Use after command 5)               */
/*                                                                                                             */
#include <stdlib.h>
#include <stdio.h>
#include "parameters.h"
#include "eval.h"
#include "gpt.h"
#include "multiMatch.h"
#include "singleMatch.h"
#include "utility.h"

#ifndef TESTMAIN
int main(int argc, char *argv[]) {
	MultiMatchContext *mmContext;
	SingleContext     *sContext;
	GptContext        *gptContext;
	FILE              *errList;

	int  command;                /* Command no */
	int  kChk = 0;                   /* k for Checked */
	char workDir[64] = WORKDIR;  /* Set work directory */

	/* For single match */
	char errListFile[MAX_FILENAME];  /* File name of output image           */
	char str[256];

	/* Image for registered and input */
	int tsCatN,                         tsNinCat;
	int rgCatN,        rgDilerN,        rgNinCat;
	int rgCatNCorrect, rgDilerNCorrect, rgNinCatCorrect;
	int dilerN;
	int rgCatNL[10],   rgNinCatL[10],  order;
#ifdef DEBUGWEIGHT
	int ix, iy; /* For debug */
#endif
	/* Check the number of arguments */
	if (argc < 2) 	{
		printf("Number of arguments is to small");
		exit(0);
	}
	command = atoi(argv[1]);

	/* Start process according to command no. */
	switch (command) {
	case 0:
		if (argc != 3) {
			printf("Number of arguments is incorrect");
			exit(0);
		}
		mmContext =	multiMatchInit(argv[2]);
		multiMatch(0, 9, 0, -1, mmContext);
		break;
	case 1:  /* Match for plural data for single registered image*/
		if (argc == 5) {
			mmContext =	multiMatchInit(argv[2]);
			multiMatch(atoi(argv[3]), atoi(argv[4]),  0, -1, mmContext);
		} else	if (argc == 7) {
			mmContext =	multiMatchInit(argv[2]);
			multiMatch(atoi(argv[3]), atoi(argv[4]),  atoi(argv[5]), atoi(argv[6]), mmContext);
		} else {
			printf("Number of arguments is incorrect");
			exit(0);
		}
		mmContext =	multiMatchInit(argv[2]);
		multiMatch(atoi(argv[3]), atoi(argv[4]),  atoi(argv[5]), atoi(argv[6]), mmContext);
		break;
	case 5: /* Make image information file */
		if (argc != 3) {
			printf("Number of arguments is incorrect");
			exit(0);
		}
		mmContext =	multiMatchInit(argv[2]);
		mkImgInfFile(mmContext);
		break;
	case 10: /* Output kNN results (Er) */
	case 11: /* Output top MAXK evaluations of results (Lt) */
	case 12: /* Output averaged evaluations (Le) */
	case 14: /* Output serials if the top of all is not correct */
	case 15: /* Output serials if kNN is not correct */
		mmContext =	multiMatchInit(argv[2]);
		if (command == 15) kChk = atoi(argv[3]);
		multiKNN(command, mmContext, kChk);
		break;
	case 20: /* Match for one to one match */
		sContext   = singleMatchInit(workDir);
		gptContext = sContext->gptContext;
		/* Registered and test image file pathes */
		sContext->tsCatN   = 1;  sContext->tsNinCat = 368;
		sContext->rgCatN   = 1;  sContext->rgNinCat = 5230;
		/* Type and file name flag */
		sContext->type     = -1;  sContext->outFlag  = "";
		switch(sContext->type) {
		case -1:
			sprintf(sContext->rgFile,     "%s/tmpImg/tdg8_0055_gray-PPT0.01.pgm", gptContext->workDir);
			sprintf(sContext->tsFile,     "%s/tmpImg/tdg8_0055_gray.pgm", gptContext->workDir);
			break;
		case 0:
			sprintf(sContext->rgFile, "%s/%s/ldg%d_%04d_gray.pgm", workDir, PGMDIR, sContext->rgCatN, sContext->rgNinCat);
			sprintf(sContext->tsFile, "%s/%s/tdg%d_%04d_gray.pgm", workDir, PGMDIR, sContext->tsCatN, sContext->tsNinCat);
			break;
		case 1:
			sprintf(sContext->rgFile, "%s%s_rg_imgInf_%d", workDir, argv[2], sContext->rgCatN);
			sprintf(sContext->tsFile, "%s%s_ts_imgInf_%d", workDir, argv[2], sContext->tsCatN);
			break;
		}
		printf("rgFile: %s\n", sContext->rgFile); printf("tsFile: %s\n", sContext->tsFile);

		/* Match */
		singleMatch(sContext);

		switch(sContext->type) {
		case -1:
			sprintf(sContext->outImgFile, "%s/tmpImg/tdg8_0055_gray-PPT0.01-FGPT.pgm", gptContext->workDir);
			save_image_file(sContext->outImgFile, sContext->rgTransImg, sContext->gptContext->nx, sContext->gptContext->ny);
			sprintf(sContext->outImgFile, "%s/tmpImg/tdg8_0055_gray-FGPT.pgm", gptContext->workDir);
			save_image_file(sContext->outImgFile, sContext->tsTransImg, sContext->gptContext->nx, sContext->gptContext->ny);
			sprintf(sContext->outImgFile, "%s/outImg/outImg_%d.pgm", gptContext->workDir, gptContext->gptDilerN);
			break;
		case 0:
			sprintf(sContext->outImgFile, "%s/outImg/out_%d_%04d_%s_%d_%04d_%1d.pgm",
					gptContext->workDir, sContext->tsCatN, sContext->tsNinCat, sContext->outFlag, sContext->rgCatN, sContext->rgNinCat, gptContext->gptDilerN);
			break;
		case 1:
			sprintf(sContext->outImgFile, "%s/outImg/out_%d_%04d_%s_%d_%04d_%1d.pgm",
					gptContext->workDir, sContext->tsCatN, sContext->tsNinCat, sContext->outFlag, sContext->rgCatN, sContext->rgNinCat, gptContext->gptDilerN);
			break;
		}
		save4Image(sContext); /* rgOrg  tsTran  tsOrg  rgTran */
#ifdef DEBUGDIFFER
			sprintf(sContext->outImgFile, "%s/outImg/out_%d_%04d_%s_%d_%04d_%1dTs.pgm",
					gptContext->workDir, sContext->tsCatN, sContext->tsNinCat, sContext->outFlag, sContext->rgCatN, sContext->rgNinCat, gptContext->gptDilerN);
			save_image_file(sContext->outImgFile, sContext->diffErImgTs, gptContext->nx, gptContext->ny);
			sprintf(sContext->outImgFile, "%s/outImg/out_%d_%04d_%s_%d_%04d_%1dTr.pgm",
					gptContext->workDir, sContext->tsCatN, sContext->tsNinCat, sContext->outFlag, sContext->rgCatN, sContext->rgNinCat, gptContext->gptDilerN);
			save_image_file(sContext->outImgFile, sContext->diffErImgTr, gptContext->nx, gptContext->ny);
#endif
#ifdef DEBUGWEIGHT
			printf("Weight of test image \n");
			for(iy = 0 ; iy < gptContext->ny ; ++iy) {
				for(ix = 0 ; ix < gptContext->nx ; ++ix) {
					printf("%5.1f ", sContext->tsImgInf->weight[ix + gptContext->nx * iy]);
					sContext->weightInt[ix + gptContext->nx * iy] = (int) sContext->tsImgInf->weight[ix + gptContext->nx * iy];
				}
				printf("\n");
			}
			printf("Direction of test image \n");
			for(iy = 0 ; iy < gptContext->ny ; ++iy) {
				for(ix = 0 ; ix < gptContext->nx ; ++ix) {
					if (sContext->tsImgInf->ang[ix + gptContext->nx * iy] == 8 ) printf("      ");
					else printf("%5d ", sContext->tsImgInf->ang[ix + gptContext->nx * iy]);
				}
				printf("\n");
			}
			sprintf(sContext->outImgFile, "%s/outImg/out_%d_%04d_%s_%d_%04d_%1dW.pgm",
					gptContext->workDir, sContext->tsCatN, sContext->tsNinCat, sContext->outFlag, sContext->rgCatN, sContext->rgNinCat, gptContext->gptDilerN);
			save_image_file(sContext->outImgFile, sContext->weightInt, gptContext->nx, gptContext->ny);
#endif
			break;
	case 21: /* GPTv1 21 execName Results/sgptError.txt : output images of error list of NN */
		sContext = singleMatchInit(workDir);
		gptContext = sContext->gptContext;
		if (argc != 4) {
			printf("Input error list\n");
			exit(0);
		}
		sprintf(sContext->execName, "%s", argv[2]);
		sprintf(errListFile, "/home/mpi/GPT%s/%s", VERSION, argv[3]);
		errList = fopen(errListFile, "r");
		if (errList == NULL) {
			printf("Error file cannot be read %s\n", errListFile);
			exit(0);
		}
		while(fgets(str, 128, errList) != NULL) {
			sscanf(str, "%d %d %d %d %d %d %d %d", &tsCatN, &tsNinCat,
					&rgCatN,  &rgNinCat, &rgDilerN, &rgCatNCorrect, &rgNinCatCorrect, &rgDilerNCorrect);
			printf("tsCatN, tsNinCat, rgCatN, rgNinCat, rgDilerN, rgCatNCorrect, rgNinCatCorrect, rgDilerNCorrect = %d, %d, %d, %d, %d, %d, %d, %d \n",
					tsCatN, tsNinCat, rgCatN, rgDilerN, rgNinCat, rgCatNCorrect, rgDilerNCorrect, rgNinCatCorrect);

			/* Best in the all categories other than correct */
			sContext->tsCatN   = tsCatN;   sContext->tsNinCat = tsNinCat;
			sContext->rgCatN   = rgCatN;   sContext->rgNinCat = rgNinCat;
			sContext->outFlag  = "D";
			sContext->type     = 0;
			switch(sContext->type) {
			case 0:
				sprintf(sContext->rgFile, "%s/%s/ldg%d_%04d_gray.pgm", workDir, PGMDIR, sContext->rgCatN, sContext->rgNinCat);
				sprintf(sContext->tsFile, "%s/%s/tdg%d_%04d_gray.pgm", workDir, PGMDIR, sContext->tsCatN, sContext->tsNinCat);
				break;
			case 1:
				sprintf(sContext->rgFile, "%s%s_rg_imgInf_%d", workDir, argv[2], sContext->rgCatN);
				sprintf(sContext->tsFile, "%s%s_ts_imgInf_%d", workDir, argv[2], sContext->tsCatN);
				break;
			}

			singleMatch(sContext);

			sprintf(sContext->outImgFile, "%s/outImg/out_%d_%04d_%s_%d_%04d_%1d.pgm",
					gptContext->workDir, sContext->tsCatN, sContext->tsNinCat, sContext->outFlag, sContext->rgCatN, sContext->rgNinCat, gptContext->gptDilerN);
			save4Image(sContext); /* rgOrg tsTrans tsOrg rgTrans */

			/* Best in the correct category */
			sContext->rgCatN   = rgCatNCorrect;  sContext->rgNinCat = rgNinCatCorrect;
			sContext->outFlag  = "S";
			switch(sContext->type) {
			case 0:
				sprintf(sContext->rgFile, "%s/%s/ldg%d_%04d_gray.pgm", workDir, PGMDIR, sContext->rgCatN, sContext->rgNinCat);
				break;
			case 1:
				sprintf(sContext->rgFile, "%s%s_rg_imgInf_%d", workDir, argv[2], sContext->rgCatN);
				break;
			}
			singleMatch(sContext);

			sprintf(sContext->outImgFile, "%s/outImg/out_%d_%04d_%s_%d_%04d_%1d.pgm",
					gptContext->workDir, sContext->tsCatN, sContext->tsNinCat, sContext->outFlag, sContext->rgCatN, sContext->rgNinCat, gptContext->gptDilerN);
			save4Image(sContext); /* rgTran  tsOrg  rgOrg  tsTran*/
		}
		break;
	case 22: /* GPTv1 22 execName Results/kNNErrorList.txt : output images of error list of kNN */
		sContext = singleMatchInit(workDir);
		gptContext = sContext->gptContext;
		if (argc != 4) {
			printf("Input error list\n");
			exit(0);
		}
		sprintf(sContext->execName, "%s", argv[2]);
		sprintf(errListFile, "/home/mpi/OpenGPT%s/%s", VERSION, argv[3]);
		errList = fopen(errListFile, "r");
		if (errList == NULL) {
			printf("Error file cannot be read %s\n", errListFile);
			exit(0);
		}
		while(fgets(str, 256, errList) != NULL) {
			sscanf(str, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d ", &tsCatN, &tsNinCat, &rgCatN,
					&(rgCatNL[0]), &(rgNinCatL[0]), &(rgCatNL[1]), &(rgNinCatL[1]), &(rgCatNL[2]), &(rgNinCatL[2]), &(rgCatNL[3]), &(rgNinCatL[3]),
					&(rgCatNL[4]), &(rgNinCatL[4]), &(rgCatNL[5]), &(rgNinCatL[5]), &(rgCatNL[6]), &(rgNinCatL[6]), &(rgCatNL[7]), &(rgNinCatL[7]),
					&(rgCatNL[8]), &(rgNinCatL[8]), &(rgCatNL[9]), &(rgNinCatL[9]));
			printf("tsCatN, tsNinCat, recogCatN = %d, %d, %d, ", tsCatN, tsNinCat, rgCatN);

			/* Best in the all categories other than correct */
			sContext->tsCatN   = tsCatN;   sContext->tsNinCat = tsNinCat;
			sContext->type     = 0;
			for (order = 0 ; order < MAXK ; ++order) {
				sContext->rgCatN   = rgCatNL[order];   sContext->rgNinCat = rgNinCatL[order];
				switch(sContext->type) {
				case 0:
					sprintf(sContext->rgFile, "%s/%s/ldg%d_%04d_gray.pgm", workDir, PGMDIR, sContext->rgCatN, sContext->rgNinCat);
					sprintf(sContext->tsFile, "%s/%s/tdg%d_%04d_gray.pgm", workDir, PGMDIR, sContext->tsCatN, sContext->tsNinCat);
					break;
				case 1:
					sprintf(sContext->rgFile, "%s%s_rg_imgInf_%d", workDir, argv[2], sContext->rgCatN);
					sprintf(sContext->tsFile, "%s%s_ts_imgInf_%d", workDir, argv[2], sContext->tsCatN);
					break;
				}

				singleMatch(sContext);

				sprintf(sContext->outImgFile, "%s/outImg/out_%1d_%04d_%1d_%1d_%1d_%04d_%1d.pgm",
						gptContext->workDir, sContext->tsCatN, sContext->tsNinCat, order, rgCatN, sContext->rgCatN, sContext->rgNinCat, gptContext->gptDilerN);
				save4Image(sContext); /* rgTran  tsOrg  rgOrg  tsTran*/
			}
		}
		break;
	case 23: /* Output dilated images */
		sContext   = singleMatchInit(workDir);
		gptContext = sContext->gptContext;
		sContext->tsCatN   = 0;  sContext->tsNinCat = 808;
		sContext->rgCatN   = 8;  sContext->rgNinCat = 4589;
		sContext->type     = 0;  sContext->outFlag  = "";
		sprintf(sContext->rgFile, "%s/%s/ldg%d_%04d_gray.pgm", workDir, PGMDIR, sContext->rgCatN, sContext->rgNinCat);
		sprintf(sContext->tsFile, "%s/%s/tdg%d_%04d_gray.pgm", workDir, PGMDIR, sContext->tsCatN, sContext->tsNinCat);
		singleMatch(sContext);

		for (dilerN = 0 ; dilerN < gptContext->nDiler ; ++dilerN) {
			sprintf(sContext->outImgFile, "%s/outImg/diler_%d_%04d_%1d.pgm", 	sContext->gptContext->workDir, sContext->rgCatN, sContext->rgNinCat, dilerN);
			printf("%s \n", sContext->outImgFile);
			if (dilerN == 0) {
				save_image_file(sContext->outImgFile, gptContext->inpImgInf->img, gptContext->nx, gptContext->ny);
			} else {
				save_image_file(sContext->outImgFile, &(gptContext->inpImgInf->imgDiler[(dilerN - 1) * gptContext->nxy]), gptContext->nx, gptContext->ny);
			}
		}
		break;
	case 24: /* Output smoothed image images */
		sContext   = singleMatchInit(workDir);
		gptContext = sContext->gptContext;
		sContext->tsCatN   = 8;  sContext->tsNinCat = 55;
		sContext->rgCatN   = 8;  sContext->rgNinCat = 4589;
		sContext->type     = 0;  sContext->outFlag  = "";
		sprintf(sContext->rgFile, "%s/%s/ldg%d_%04d_gray.pgm", workDir, PGMDIR, sContext->rgCatN, sContext->rgNinCat);
		sprintf(sContext->tsFile, "%s/%s/tdg%d_%04d_gray.pgm", workDir, PGMDIR, sContext->tsCatN, sContext->tsNinCat);
		/* We do not use matching. Only use arrays */
		singleMatch(sContext);
		dilerN = 2; /* The number of Smoothing */

		imgInfSmoothD(sContext->tsImgInf->img, sContext->tsImgInf->can, dilerN, sContext->rgImgInf->can, gptContext->nx, gptContext->ny);
		int ixy;
		for (ixy = 0 ; ixy < gptContext->nxy ; ++ixy) sContext->tsImgInf->img[ixy] = (int) sContext->tsImgInf->can[ixy];

		sprintf(sContext->outImgFile, "%s/outImg/smooth_%d_%04d_%1d.pgm", workDir, sContext->tsCatN, sContext->tsNinCat, dilerN);
		printf("%s \n", sContext->outImgFile);
		save_image_file(sContext->outImgFile, sContext->tsImgInf->img, gptContext->nx, gptContext->ny);
		break;
	}
	return 0;
}

#else
#include "utility.h"
	/*
	 * Usage : ./gptMain type coef nIter
	 *        ex.:  ./gptMain 0 0.3 10
	 */
	int main(int argc, char *argv[]) {
		int type = atoi(argv[1]), nIter = atoi(argv[3]);
		double coef; sscanf(argv[2], "%lf", &coef);
		int nx = COL, ny = ROW;
		int pos;
		char workDir[64] = WORKDIR;  /* Set work directory */
		/* For single match */
		char rgFile[MAX_FILENAME];      /* File name of registered images      */
		char outImgFile[MAX_FILENAME];  /* File name of output image           */

		int     *img      = (int *) malloc(sizeof(int) * nx * ny);
		double  *inImgD   = (double *) malloc(sizeof(double) * nx * ny);
		double  *outImgD  = (double *) malloc(sizeof(double) * nx * ny);
		double  *workImgD = (double *) malloc(sizeof(double) * nx * ny);

		sprintf(rgFile,     "%s/dilImage/input.pgm", workDir);
		//	sprintf(rgFile,     "%s/dilImage/pulse.pgm", workDir);
		sprintf(outImgFile, "%s/dilImage/out_%d_%f_%d.pgm", workDir, type, coef, nIter);
		load_image_file(rgFile, img, nx, ny, 0);
		for (pos = 0 ; pos < nx * ny ; ++pos) {
			inImgD[pos] = img[pos];
		}

		dilation(type, inImgD, outImgD, workImgD, nx, ny, coef, nIter);
		for (pos = 0 ; pos < nx * ny ; ++pos) {
			img[pos] = (int) outImgD[pos];
		}
		save_image_file(outImgFile, img, nx, ny);
		return 0;
	}
#endif
