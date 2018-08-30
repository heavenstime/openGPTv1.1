/*
 * gpt.c
 *
 *  Created on: 2012/11/02
 *      Author: yamasita
 */

/* Parameters in program                               */
/*   gammaList : List of parameters of Gauss integral  */
/*                                                                   */
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "parameters.h"
#include "eval.h"
#include "gpt.h"
#include "utility.h"

#define C11 0
#define C21 1
#define C31 2
#define C12 3
#define C22 4
#define C32 5
#define C13 6
#define C23 7
#define C33 8

void gptInvTransformImage(double *gpt, int *inImg, int *outImg, GptContext *gptContext);
void gptInvTransformImageD(double *gpt, double *inImg, double *outImg, GptContext *gptContext);
void gptInit(double *gpt);
void gptCopy(double *inGpt, double *outGpt);
void gptTransformPoint(double *gpt, double *inP, double *outP);
void gptTransformGpt(double *gpt1, double *gpt2, double *outGpt);
int  gptInverse(double *gpt, double *iGpt);

//#define DEBUG
//#define DEBUGITER
//#define DEBUGPrG
//#define DEBUGEVAL

//#define NODEFCAN

/* Allocate region for calculation of GPT */
void gptCorInit(int nx, int ny, int nExt, int nGamma, double *gammaList, GptContext *gptContext) {
	int nxy, nxyDirNdGamma, nxy2;
	int lGamma;

	/* Parameter independent of content of image */
	gptContext->nx          = nx;
	gptContext->ny          = ny;
	gptContext->nExt        = nExt;;
	gptContext->cx          = nx / 2;
	gptContext->cy          = ny / 2;
	gptContext->nGamma      = nGamma;    /* # of element in gammaList */
	gptContext->gammaList   = (double *) malloc(sizeof(double) * nGamma);
	gptContext->nDir        = DIRECTION;
	gptContext->lTrunc      = 2.5;  /* Coefficient to parameter of Gauss integral for its truncation */
	gptContext->nDiffEr     = NDIFFER;
	gptContext->inverseChk  = INVERSECHK;
	gptContext->transChk    = TRANSCHK;
	gptContext->smooth      = INPUTSMOOTH;

	/* Constant decided from parameters */
	gptContext->nxy           = nxy  = nx * ny;
	gptContext->nxy2          = nxy2 = (nx + 2) * (ny + 2);
	gptContext->nxyDir        = nxy * gptContext->nDir;
#ifdef DIRMATCH
	gptContext->nDirNd        = DIRECTION + 1;
#else
	gptContext->nDirNd        = 1;
#endif
	gptContext->nxyDirNd      = nxy * gptContext->nDirNd;
	gptContext->nxyDirNdGamma = nxyDirNdGamma = gptContext->nxyDirNd * gptContext->nGamma;
	gptContext->nst           = (2 * nx - 1) * (2 * ny - 1);
	gptContext->dDiagSq       = nx * nx + ny * ny;
	gptContext->dDiag         = sqrt((double) gptContext->dDiagSq);
	/* Tables which do not depend on content of image */
	gptContext->xFull       = (GptTblLine *) malloc(sizeof(GptTblLine) * gptContext->nGamma);
	gptContext->yFull       = (GptTblLine *) malloc(sizeof(GptTblLine) * gptContext->nGamma);
	gptContext->xTrunc      = (GptTblLine *) malloc(sizeof(GptTblLine) * gptContext->nGamma);
	gptContext->yTrunc      = (GptTblLine *) malloc(sizeof(GptTblLine) * gptContext->nGamma);
	for (lGamma = 0 ; lGamma < gptContext->nGamma ; ++lGamma) {
		gptContext->gammaList[lGamma]   = gammaList[lGamma];

		gptContext->xFull[lGamma].gauss = (double *) malloc(sizeof(double) * (nx + nx - 1));
		gptContext->xFull[lGamma].st    = (int *) malloc(sizeof(int) * nx);
		gptContext->xFull[lGamma].gSt   = (int *) malloc(sizeof(int) * nx);
		gptContext->xFull[lGamma].gEnd  = (int *) malloc(sizeof(int) * nx);

		gptContext->yFull[lGamma].gauss = (double *) malloc(sizeof(double) * (ny + ny - 1));
		gptContext->yFull[lGamma].st    = (int *) malloc(sizeof(int) * ny);
		gptContext->yFull[lGamma].gSt   = (int *) malloc(sizeof(int) * ny);
		gptContext->yFull[lGamma].gEnd  = (int *) malloc(sizeof(int) * ny);

		gptContext->xTrunc[lGamma].gauss = (double *) malloc(sizeof(double) * (nx + nx - 1));
		gptContext->xTrunc[lGamma].st    = (int *) malloc(sizeof(int) * nx);
		gptContext->xTrunc[lGamma].gSt   = (int *) malloc(sizeof(int) * nx);
		gptContext->xTrunc[lGamma].gEnd  = (int *) malloc(sizeof(int) * nx);

		gptContext->yTrunc[lGamma].gauss = (double *) malloc(sizeof(double) * (ny + ny - 1));
		gptContext->yTrunc[lGamma].st    = (int *) malloc(sizeof(int) * ny);
		gptContext->yTrunc[lGamma].gSt   = (int *) malloc(sizeof(int) * ny);
		gptContext->yTrunc[lGamma].gEnd  = (int *) malloc(sizeof(int) * ny);
	}
	/* Tables for original algorithm */
	gptContext->gDist     = (double *) malloc(sizeof(double) * nxy);
	gptContext->gwt       = (double *) malloc(sizeof(double) * nxy);
	/* Work Region */
	gptContext->gh0       = (double *) malloc(sizeof(double) * nxyDirNdGamma);
	gptContext->gh1       = (double *) malloc(sizeof(double) * nxyDirNdGamma);
	gptContext->gh2       = (double *) malloc(sizeof(double) * nxyDirNdGamma);
	gptContext->atrImg    = (int *) malloc(sizeof(int) * nxy);
	gptContext->atrAng    = (int *) malloc(sizeof(int) * nxy);
	gptContext->atrCan    = (double *) malloc(sizeof(double) * nxy);
	gptContext->atrWeight = (double *) malloc(sizeof(double) * nxy);
	gptContext->invAtrImg = (int *) malloc(sizeof(int) * nxy);
	gptContext->invAtrCan = (double *) malloc(sizeof(double) * nxy);
	gptContext->atrImgNew = (int *) malloc(sizeof(int) * nxy);
	gptContext->workImgD1 = (double *) malloc(sizeof(double) * nxy2);
	gptContext->workImgD2 = (double *) malloc(sizeof(double) * nxy2);
	gptContext->workImgD3 = (double *) malloc(sizeof(double) * nxy2);
	gptContext->workImgD4 = (double *) malloc(sizeof(double) * nxy2);
	nndegdInit(nx, ny, gptContext->nndegdContext);
#ifdef DEBUGDIFFER
	gptContext->diffErImg = (int *) malloc(sizeof(int) * gptContext->nxy); 	/* To keep eroded difference image */
#endif
	gptContext->debug = 0;
}
/* Calculation of GPT correlation */
int calGptCor(GptContext *gptContext) {
	int     *atrImg    = gptContext->atrImg;     /* Affine transformed image  */
	int     *atrAng    = gptContext->atrAng;
	double  *atrCan    = gptContext->atrCan;
	double  *atrWeight = gptContext->atrWeight;
	int     *atrImgNew = gptContext->atrImgNew ; /* Newly affine transformed image */
	int     *invAtrImg = gptContext->invAtrImg;
	double  *invAtrCan = gptContext->invAtrCan;
	ImgInf  *inpImgInf = gptContext->inpImgInf;
	ImgInf  *tgtImgInf = gptContext->tgtImgInf;

	NndegdContext *nndegdContext = gptContext->nndegdContext;

	double *gpt    = gptContext->gpt;
	double *iGpt   = gptContext->iGpt;
	double *gptNew = gptContext->gptNew;

	int     nx      = gptContext->nx;
	int     ny      = gptContext->ny;
	int     nxy     = gptContext->nxy;
	int     nxyDir  = gptContext->nxyDir;

	double  corNew, selfCor; /* Correlations */
	double  smpEval, gptEval;
	int     dilerN = 0;
	int     iter = 0;        /* # of iterations */
	int     transType = 0;
	int     gptStatus = 0;
	int     corStatus, invStatus;

	smpEval = gptEval = imgCor(gptContext->inpImgInf->can, gptContext->tgtImgInf->can, nxy);
#ifdef DEBUGITER
	if (gptContext->debug) {
		printf("iter = %d  cor = %f  \n", iter, gptEval);
	}
#endif

#ifndef NOGPT
	/* Initialize Affine transform and data */
	gptInit(gpt);
	imgCopy(inpImgInf->img, atrImg, nxy);
	imgCopy(inpImgInf->ang, atrAng, nxy);
	canCopy(inpImgInf->can, atrCan, nxy);
	gptCopy(gpt, gptNew);
	for ( ; iter < MAX_ITER ; ++iter) {

		/* Obtain the parameter of Gauss integral */
		gptContext->nndegd = nndegd(NNDEGDTYPE, atrAng, tgtImgInf->ang, nndegdContext, tgtImgInf->nndegdTbl);
		if (gptContext->nndegd == NNDEGDFAIL) {
			gptStatus = GPTFAIL;
			break;
		}

		/* Obtain Affine transform for GPT correlation */
		switch(gptContext->gptType) {
		case ALLGAT:
			transType = GPTGAT;
			break;
		case ALTGATPPT:
			transType = (iter % 2 == 0) ? GPTGAT : GPTPPT;
			break;
		case GATPPTAFTGAT:
			transType = (iter < TIMESFSTGAT || (iter % 2 != TIMESFSTGAT % 2)) ? GPTGAT : GPTPPT;
			break;
		}

		corStatus = gptCor(transType, atrAng, atrCan, tgtImgInf->ang, tgtImgInf->can, gptContext->nndegd, gptNew, gptContext);
#ifdef DEBUGITER
		if (gptContext->debug == 1) gptPr(gptNew, "New");
#endif
		if (corStatus < 0) { /* Error in calculation of GPT */
			gptStatus = GPTFAIL;
			break;
		}
		/* Affine transform image by the above Affine transform */
		invStatus = gptInverse(gptNew, iGpt);
		if (invStatus != 0) {
			gptStatus = GPTFAIL;
			break;
		}

		gptInvTransformImage(iGpt, inpImgInf->img, atrImgNew, gptContext);  /* Affine transform of image */
		/* imgPr(inpImgInf->img, atrImgNew, nx, ny); */
		/* save_image_file("gptImg.pgm", atrImgNew, nx, ny); */

		/* Obtain canonical image of newly transformed image */
		/* and calculation of correlation */
		defcan(atrImgNew, atrCan, nxy);
		corNew = imgCor(atrCan, tgtImgInf->can, nxy);

#ifdef DEBUGITER
		if (gptContext->debug) {
			printf("iter = %d  cor = %f nndegd = %f gamma = %lf \n", iter + 1, corNew, gptContext->nndegd,
					0.5 / ( WGT * WGT * gptContext->nndegd * gptContext->nndegd));
		}
#else
		/* If correlation is not decreased, exit the loop */
		if (corNew < STOP_ITER * gptEval) break;
#endif
		gptEval = corNew;
		imgCopy(atrImgNew, atrImg, nxy); /* Since the new one is better that old one, the new is copied to artImg */
		gptCopy(gptNew, gpt);

		/* Calculate of edge angle image from a newly transformed image */
		roberts8(atrImg, atrAng, gptContext->nx, gptContext->ny);
	}
#endif

	/* Output final evaluation and iteration */
	/* Set initial values without dilation/erosion */
#ifndef NOGPT
	gptContext->smpDilerN = gptContext->gptDilerN = 0;
	switch(gptContext->evalType) {
	case EVALCOR:
		gptContext->smpEval = smpEval;
		gptContext->gptEval = gptEval;
		break;
	case EVALSSD:
		break;
	case EVALNNDEGDWEIGHT:
	case EVALNNDEGDSQRWEIGHT:
		calWeight(atrImg, atrWeight, nx, ny, gptContext->workImgD1, gptContext->workImgD2, gptContext->workImgD3, gptContext->workImgD4);
	case EVALNNDEGD:
	case EVALNNDEGDVAR:
	case EVALNNDEGDSQR:
	case EVALNNDEGDSQRT:
		gptContext->smpEval = - nndegdEv(gptContext->evalDnnType, gptContext->evalType, inpImgInf->ang, tgtImgInf->ang, inpImgInf->weight, tgtImgInf->weight, nndegdContext, tgtImgInf->nndegdTbl);
		gptContext->gptEval = - nndegdEv(gptContext->evalDnnType, gptContext->evalType, atrAng,         tgtImgInf->ang, atrWeight,         tgtImgInf->weight, nndegdContext, tgtImgInf->nndegdTbl);
		// gptContext->gptEval = - nndegd(gptContext->evalDnnType, atrAng,         tgtImgInf->ang, nndegdContext, tgtImgInf->nndegdTbl);
		break;
	case EVALDIFFER:
		gptContext->smpEval = - diffImgEr(inpImgInf->img, tgtImgInf->img, gptContext->nDiffEr, nx, ny, gptContext->workImgD1, gptContext->workImgD2, gptContext->workImgD3, gptContext->workImgD4, gptContext->diffErImg);
		gptContext->gptEval = - diffImgEr(        atrImg, tgtImgInf->img, gptContext->nDiffEr, nx, ny, gptContext->workImgD1, gptContext->workImgD2, gptContext->workImgD3, gptContext->workImgD4, gptContext->diffErImg);
		break;
	case EVALWDCH:
		gptContext->smpEval = wdchEv(inpImgInf->ang, tgtImgInf->ang, nx, ny);
		gptContext->gptEval = wdchEv(        atrAng, tgtImgInf->ang, nx, ny);
		break;
	}
#ifdef DEBUGEVAL
	printf("dilerN = %d, smpEval = %f, gptEval = %f \n", dilerN, gptContext->smpEval, gptContext->gptEval);
#endif
	if (gptContext->gptEval > DILCUTTHRESH) { /* Dileration is not done if correlation is small */
		for (dilerN = 1 ; dilerN < gptContext->nDiler ; ++dilerN) {
			/* Output value for recognition */
			switch(gptContext->evalType) {
			case EVALCOR:
				smpEval = imgCor(inpImgInf->can, &(tgtImgInf->canDiler[(dilerN - 1) * nxy]), nxy);
				defcan(atrImg, atrCan, nxy);
				gptEval = imgCor(atrCan,         &(tgtImgInf->canDiler[(dilerN - 1) * nxy]), nxy);
				break;
			case EVALNNDEGDWEIGHT:
			case EVALNNDEGDSQRWEIGHT:
				calWeight(atrImg, atrWeight, nx, ny, gptContext->workImgD1, gptContext->workImgD2, gptContext->workImgD3, gptContext->workImgD4);
			case EVALNNDEGD:
			case EVALNNDEGDVAR:
			case EVALNNDEGDSQR:
			case EVALNNDEGDSQRT:
				smpEval = - nndegdEv(gptContext->evalDnnType, gptContext->evalType, inpImgInf->ang, &(tgtImgInf->angDiler[(dilerN - 1) * nxy]), inpImgInf->weight, tgtImgInf->weight,                          nndegdContext, &(tgtImgInf->nndegdTblDiler[(dilerN - 1) * nxyDir]));
				gptEval = - nndegdEv(gptContext->evalDnnType, gptContext->evalType, atrAng,         &(tgtImgInf->angDiler[(dilerN - 1) * nxy]), atrWeight,         &(tgtImgInf->weightDiler[(dilerN - 1) * nxy]), nndegdContext, &(tgtImgInf->nndegdTblDiler[(dilerN - 1) * nxyDir]));
				break;
			}
			if (smpEval > gptContext->smpEval) {
				gptContext->smpEval   = smpEval;
				gptContext->smpDilerN = dilerN;
			}
			if (gptEval > gptContext->gptEval) {
				gptContext->gptEval   = gptEval;
				gptContext->gptDilerN = dilerN;
			}
		}
#ifdef DEBUGEVAL
		printf("dilerN = %d, smpEval = %f, gptEval = %f \n", dilerN, smpEval, gptEval);
#endif
	}
#else
	gptContext->gptEval = gptEval;
	switch(gptContext->evalType) {
	case EVALCOR:
		gptContext->smpEval = smpEval;
		break;
	case EVALSSD:
		break;
	case EVALNNDEGDWEIGHT:
	case EVALNNDEGDSQRWEIGHT:
	case EVALNNDEGD:
	case EVALNNDEGDVAR:
	case EVALNNDEGDSQR:
	case EVALNNDEGDSQRT:
		gptContext->smpEval = - nndegdEv(gptContext->evalDnnType, gptContext->evalType, inpImgInf->ang, tgtImgInf->ang, inpImgInf->weight, tgtImgInf->weight, nndegdContext, tgtImgInf->nndegdTbl);
		break;
	case EVALDIFFER:
		gptContext->smpEval = - diffImgEr(inpImgInf->img, tgtImgInf->img, gptContext->nDiffEr, nx, ny, gptContext->workImgD1, gptContext->workImgD2, gptContext->workImgD3, gptContext->workImgD4, gptContext->diffErImg);
		break;
	case EVALWDCH:
		gptContext->smpEval = wdchEv(inpImgInf->ang, tgtImgInf->ang, nx, ny);
		break;
	}
#ifdef DEBUGEVAL
	printf("dilerN = %d, smpEval = %f, gptEval = %f \n", dilerN, gptContext->smpEval, gptContext->gptEval);
#endif
	if (gptContext->gptEval > DILCUTTHRESH) { /* Dileration is not done if correlation is small */
		for (dilerN = 1 ; dilerN < gptContext->nDiler ; ++dilerN) {
			/* Output value for recognition */
			switch(gptContext->evalType) {
			case EVALCOR:
				smpEval = imgCor(inpImgInf->can, &(tgtImgInf->canDiler[(dilerN - 1) * nxy]), nxy);
				break;
			case EVALNNDEGDWEIGHT:
			case EVALNNDEGDSQRWEIGHT:
			case EVALNNDEGD:
			case EVALNNDEGDVAR:
			case EVALNNDEGDSQR:
			case EVALNNDEGDSQRT:
				smpEval = - nndegdEv(gptContext->evalDnnType, gptContext->evalType, inpImgInf->ang, &(tgtImgInf->angDiler[(dilerN - 1) * nxy]), inpImgInf->weight, tgtImgInf->weight,                          nndegdContext, &(tgtImgInf->nndegdTblDiler[(dilerN - 1) * nxyDir]));
				break;
			}
			if (smpEval > gptContext->smpEval) {
				gptContext->smpEval   = smpEval;
				gptContext->smpDilerN = dilerN;
			}
		}
	}
#endif

	/* Inverse transformation of transformation (to check bound) */
	if (gptContext->inverseChk == 1) {
		gptInvTransformImage(gpt, atrImg, invAtrImg, gptContext);  /* Affine transform of image */
		defcan(invAtrImg, invAtrCan, nxy);
		selfCor = imgCor(invAtrCan, inpImgInf->can, nxy);
		// printf("SelfCor = %f \n", selfCor);
		switch(gptContext->evalType) {
		case EVALCOR:
			gptContext->gptEval *= selfCor;
			break;
		default:
			if (selfCor < SELFCORLIMIT) {
				gptContext->gptEval   = -NNDEGDFAIL;
				gptContext->gptDilerN = 0;
			}
			break;
		}
	}
	/* Transform check (to check bound) */
	if (gptContext->transChk == 1) {
		double a11 = gpt[C11], a12 = gpt[C12], a21 = gpt[C21], a22 = gpt[C22];
		double tr  = a11 * a11 + a12 * a12 + a21 * a21 + a22 * a22;
		double det = (a11 * a11 + a12 * a12) * (a21 * a21 + a22 * a22) - (a11 * a21 + a12 * a22) * (a11 * a21 + a12 * a22);
		if (det < TRANSLIMIT * tr * tr * 0.25) {
			gptContext->gptEval   = -NNDEGDFAIL;
			gptContext->gptDilerN = 0;
		}
	}

	gptContext->iter = iter;
	return gptStatus;
}

/* Determination of optimal GPT components */
/* that yield the maximal correlation value */
/* gptType 0:GAT, 1:PPT, 2:(GAT + PPT) / 2 */
int gptCor(int transType, int *g_ang1, double *g_can1, int *g_ang2, double *g_can2,
		double nndegd, double *gpt, GptContext *gptContext) {
	int nx = gptContext->nx, ny = gptContext->ny, nDirNd = gptContext->nDirNd;
	int nxy = gptContext->nxy, nxyDirNd = gptContext->nxyDirNd;
	int nGamma = gptContext->nGamma;
	int x, y, x1, y1, x2, y2, pos, tPos1, tPos2, loop;
	int taPos1, taPos2, lGamma;
	double gamma, *gDist = gptContext->gDist, *gwt = gptContext->gwt;
	double wGamma1, wGamma2;
	double gOne, gx1, gy1, gx2, gy2;
	double gx1x1, gx1y1, gy1y1,	gx1x2, gx1y2, gy1x2, gy1y2;
	double pcor, t0, tx2, ty2;
	double tIn12, tIn11, gIn12x1x1, gIn12x1y1, gIn12y1y1, gIn11x1, gIn11y1, gIn12x1, gIn12y1;
	double gIn11, gIn12, v1, v2;
	double V11, V12, V21, V22, U11, U12, U21, U22, Uinv11, Uinv12, Uinv21, Uinv22, det;
	double *gH0 = gptContext->tgtImgInf->gH0, *gH1x = gptContext->tgtImgInf->gH1x, *gH1y = gptContext->tgtImgInf->gH1y;
	double *gammaList = gptContext->gammaList;
	double gcH0, gcH1x, gcH1y;
	double gptNew[9], gptGat[9], gptPpt[9];
	double dx1, dy1, dx2, dy2, dx, dy, cx = gptContext->cx, cy = gptContext->cy;

	gamma  = 0.5 / ( WGT * WGT * nndegd * nndegd);

	switch(gptContext->gptMethod) {
	case GPTORG: /* A little modified Original GPT algorithm */
		/* Make Gauss function */
		for (pos = 0 ; pos < nxy ; ++pos) {
			gwt[pos] = gamma * exp(-gamma * gDist[pos]);
		}
		/* Gaussian weigthed mean values */
		/* Initialize for GPT */
		gOne = gx1 = gy1 = gx2 = gy2 = 0.0;
		gx1x1 = gx1y1 = gy1y1 = 0.0;
		gx1x2 = gx1y2 = gy1x2 = gy1y2 = 0.0;
		/* Initialize for PPT */
		gIn12x1x1 = gIn12x1y1 = gIn12y1y1 = gIn11x1 = gIn11y1 = gIn12x1 = gIn12y1 = 0.0;
		gIn11 = gIn12 = 0.0;
		for (y1 = 0 ; y1 < ny ; ++y1) {
			dy1 = y1 - cy;
			for (x1 = 0 ; x1 < nx ; ++x1) {
				dx1 = x1 - cx;
				t0 = 0.0; tx2 = 0.0; ty2 = 0.0;
				for (y2 = 0 ; y2 < ny ; ++y2) {
					dy2 = y2 - cy;
					for (x2 = 0 ; x2 < nx ; ++x2) {
						dx2 = x2 - cx;
#ifdef DIRMATCH
						if (g_ang1[y1 * nx + x1] == g_ang2[y2 * nx + x2]) {
#else
							//	if (g_can1[y1 * nx + x1] < 0.0 && g_can2[y2 * nx + x2] < 0.0) {
							if (1) {
#endif
								pcor  = gwt[abs(y2 - y1) * nx + abs(x2 - x1)] * g_can1[y1 * nx + x1] * g_can2[y2 * nx + x2];
								t0     += pcor;
								tx2    += pcor * dx2;
								ty2    += pcor * dy2;
#ifdef DIRMATCH
							}
#else
						}
#endif
					}
				}
				gOne  += t0;
				if (transType == GPTGAT || transType == GPTGATPPT) {
					gx2 += tx2; gy2 += ty2;
					gx1   += t0 * dx1;
					gy1   += t0 * dy1;
					gx1x1 += t0 * dx1 * dx1;
					gx1y1 += t0 * dx1 * dy1;
					gy1y1 += t0 * dy1 * dy1;
					gx1x2 += tx2 * dx1;
					gx1y2 += ty2 * dx1;
					gy1x2 += tx2 * dy1;
					gy1y2 += ty2 * dy1;
				}
				if (transType == GPTPPT || transType == GPTGATPPT) {
					tIn11      = (dx1 * dx1 + dy1 * dy1) * t0;
					tIn12      =  dx1 * tx2 + dy1 * ty2;
					gIn12x1x1 += tIn12 * dx1 * dx1;
					gIn12x1y1 += tIn12 * dx1 * dy1;
					gIn12y1y1 += tIn12 * dy1 * dy1;
					gIn11x1   += tIn11 * dx1;
					gIn11y1   += tIn11 * dy1;
					gIn12x1   += tIn12 * dx1;
					gIn12y1   += tIn12 * dy1;
					gIn11     += tIn11;
					gIn12     += tIn12;

				}
			}
		}
		break;
	case GPTTMPLT:  /* Algorithm using template */
		if (gamma >= gammaList[0]) {
			/* gamma is larger than the largest in gammaList */
			tPos1 = 0; wGamma1 = 1.0;
			tPos2 = 0; wGamma2 = 0.0;
		} else {
			/* Set gamma to the smallest in gammaList */
			tPos1 = nxyDirNd * (nGamma - 1) ; wGamma1 = 1.0;
			tPos2 = nxyDirNd * (nGamma - 1) ; wGamma2 = 0.0;
			/* もう少し大きい場合を大きい方から調べる */
			for (lGamma = 1 ; lGamma < nGamma ; ++lGamma) {
				if (gamma >= gammaList[lGamma]) {
					tPos1 = nxyDirNd * (lGamma - 1);
					tPos2 = nxyDirNd * lGamma;
					wGamma1 = (gamma - gammaList[lGamma]) / (gammaList[lGamma - 1] - gammaList[lGamma]);
					wGamma2 = 1.0 - wGamma1;
					break;
				}
			}
		}
#ifdef DEBUG
		printf("gamma = %f lG = (%d, %d) wGamma = (%f, %f) \n", gamma, tPos1 / nxyDirNd, tPos2 /nxyDirNd, wGamma1, wGamma2);
#endif
		pos = 0;
		gOne  = 0.0; gx1   = 0.0; gy1   = 0.0;
		gx1x1 = 0.0; gx1y1 = 0.0; gy1y1 = 0.0;
		gx2   = 0.0; gy2   = 0.0;
		gx1x2 = 0.0; gy1x2 = 0.0; gx1y2 = 0.0; gy1y2 = 0.0;
		gIn12x1x1 = gIn12x1y1 = gIn12y1y1 = gIn11x1 = gIn11y1 = gIn12x1 = gIn12y1 = 0.0;
		gIn11 = gIn12 = 0.0;
		for (y = 0 ; y < ny ; ++y) {
			dy = y - gptContext->cy;
			for (x = 0 ; x < nx ; ++x) {
				dx = x - gptContext->cx;
#ifdef DIRMATCH
				taPos1 = tPos1 + g_ang1[pos];
				taPos2 = tPos2 + g_ang1[pos];
#else
				taPos1 = tPos1;
				taPos2 = tPos2;
#endif
				gcH0  = g_can1[pos] * (wGamma1 * gH0[taPos1]  + wGamma2 * gH0[taPos2]);
				gcH1x = g_can1[pos] * (wGamma1 * gH1x[taPos1] + wGamma2 * gH1x[taPos2]);
				gcH1y = g_can1[pos] * (wGamma1 * gH1y[taPos1] + wGamma2 * gH1y[taPos2]);
				gOne  += gcH0;
				if (transType != 1) {
					gx1   += dx * gcH0;
					gy1   += dy * gcH0;
					gx1x1 += gcH0 * dx * dx;
					gx1y1 += gcH0 * dx * dy;
					gy1y1 += gcH0 * dy * dy;
					gx2   += gcH1x;
					gy2   += gcH1y;
					gx1x2 += gcH1x * dx;
					gy1x2 += gcH1x * dy;
					gx1y2 += gcH1y * dx;
					gy1y2 += gcH1y * dy;
				}
				if (transType > 0) {
					tIn11      = (dx * dx + dy * dy) * gcH0;
					tIn12      =  dx * gcH1x + dy * gcH1y;
					gIn12x1x1 += tIn12 * dx * dx;
					gIn12x1y1 += tIn12 * dx * dy;
					gIn12y1y1 += tIn12 * dy * dy;
					gIn11x1   += tIn11 * dx;
					gIn11y1   += tIn11 * dy;
					gIn12x1   += tIn12 * dx;
					gIn12y1   += tIn12 * dy;
					gIn11     += tIn11;
					gIn12     += tIn12;
				}
				++pos; tPos1 += nDirNd; tPos2 += nDirNd;
			}
		}
		break;
	default:
		printf("Error of set of Gpt method.\n");
		return -1;
		break;
	}
	if (fabs(gOne) < EPS) {
		printf("GPT calculation failure by zero sum!!!\n");
		return -1;
	}

	/* Calc GPT for GPT */
	if (transType == GPTGAT || transType == GPTGATPPT) {
		/* bar(x_2 x_1^T) - bar(x_2) bar(x_1)^T / var(1) */
		V11 = gx1x2 - gx1 * gx2 / gOne;
		V12 = gy1x2 - gy1 * gx2 / gOne;
		V21 = gx1y2 - gx1 * gy2 / gOne;
		V22 = gy1y2 - gy1 * gy2 / gOne;

		/* bar(x_1 x_1^T) - bar(x_1) bar(x_1)^T / bar(1) */
		U11 = gx1x1 - gx1 * gx1 / gOne;
		U12 = U21 = gx1y1 - gy1 * gx1 / gOne;
		U22 = gy1y1 - gy1 * gy1 / gOne;
		det = U11 * U22 - U12 * U21;
		if (fabs(det) < EPS) {
			printf("GPT calculation failure by zero det of U !!! det = %f\n", det);
			return -1;
		}
		Uinv11 = U22 / det;
		Uinv12 = Uinv21 = -U12 / det;
		Uinv22 = U11 / det;

#ifdef DEBUGPrG
		printf("g11 : %f %f %f,  g22 : %f %f %f \n", gx1x1, gx1y1, gy1y1, gx2x2, gx2y2, gy2y2);
#endif
		gptGat[C11] = V11 * Uinv11 + V12 * Uinv21;
		gptGat[C21] = V21 * Uinv11 + V22 * Uinv21;
		gptGat[C12] = V11 * Uinv12 + V12 * Uinv22;
		gptGat[C22] = V21 * Uinv12 + V22 * Uinv22;
		gptGat[C13] = (gx2 - gptGat[C11] * gx1 - gptGat[C12] * gy1) / gOne;
		gptGat[C23] = (gy2 - gptGat[C21] * gx1 - gptGat[C22] * gy1) / gOne;
		gptGat[C31] = 0.0; gptGat[C32] = 0.0; gptGat[C33] = 1.0;
		// gptPr(gptGat, "GAT");
	}

	/* Calc GPT for PPT */
	if (transType == GPTPPT || transType == GPTGATPPT) {
		U11 =       gIn12x1x1;
		U12 = U21 = gIn12x1y1;
		U22 =       gIn12y1y1;

		/* U^{-1} */
		det = U11 * U22 - U21 * U12;
		if (fabs(det) < EPS) {
			printf("GPT calculation failure by zero det of V for Both side !!! det = %f\n", det);
			return -1;
		}
		Uinv11 =  U22 / det;
		Uinv12 = -U12 / det;
		Uinv21 = -U21 / det;
		Uinv22 =  U11 / det;

		/* V = bar(<x_1 x_1> x_1) - bar(<x_1 x_2> x_1) */
		v1 = gIn11x1 - gIn12x1;
		v2 = gIn11y1 - gIn12y1;
		/* printf("U : %f %f %f,  v : %f %f\n", U11, U12, U22, v1, v2); */

		gptPpt[C11] = 1.0; gptPpt[C12] = 0.0; gptPpt[C13] = 0.0;
		gptPpt[C21] = 0.0; gptPpt[C22] = 1.0; gptPpt[C23] = 0.0;
		gptPpt[C31] = Uinv11 * v1 + Uinv12 * v2;
		gptPpt[C32] = Uinv21 * v1 + Uinv22 * v2;
		gptPpt[C33] = 1.0;
		// gptPr(gptPpt, "PPT");
	}

	/* Update of GPT components */
	switch(transType) {
	case GPTGAT:
		gptTransformGpt(gptGat, gpt, gptNew); break;
	case GPTPPT:
		gptTransformGpt(gptPpt, gpt, gptNew); break;
	case GPTGATPPT:
		for (loop = 0 ; loop < 9 ; ++loop) gptGat[loop] = 0.5 * (gptGat[loop] + gptPpt[loop]);
		gptTransformGpt(gptGat, gpt, gptNew); break;
	}
	// gptPr(gptNew, "New");
	gptCopy(gptNew, gpt);
	return 0;
}

/* Make GPT calculation table　*/
void mkGptTbl(GptContext *gptContext) {
	int x, y, pos, lGamma;
	int nx = gptContext->nx, ny = gptContext->ny;
	double gamma;
	/* Table of square coordinates (Recently it makes slower?) */
	pos = 0;
	for (y = 0 ; y < ny ; ++y) {
		for (x = 0 ; x < nx ; ++x) {
			gptContext->gDist[pos]  = x * x + y * y;
			++pos;
		}
	}
	for (lGamma = 0 ; lGamma < gptContext->nGamma ; ++lGamma) {
		gamma = gptContext->gammaList[lGamma];
		/* x full integral */
		gptContext->xFull[lGamma].leng     = nx;
		gptContext->xFull[lGamma].intgLeng = nx - 1;
		gptContext->xFull[lGamma].gamma    = gamma;
		mkGptTblLine(& (gptContext->xFull[lGamma]) );

		/* y full integral */
		gptContext->yFull[lGamma].leng     = ny;
		gptContext->yFull[lGamma].intgLeng = ny - 1;
		gptContext->yFull[lGamma].gamma    = gamma;
		mkGptTblLine(& (gptContext->yFull[lGamma]) );

		/* x truncated integral */
		gptContext->xTrunc[lGamma].leng     = nx;
		gptContext->xTrunc[lGamma].intgLeng = calIntgLeng(nx, gamma, gptContext->lTrunc);
		gptContext->xTrunc[lGamma].gamma    = gamma;
		mkGptTblLine(& (gptContext->xTrunc[lGamma]) );

		/* y truncated integreal */
		gptContext->yTrunc[lGamma].leng     = ny;
		gptContext->yTrunc[lGamma].intgLeng = calIntgLeng(ny, gamma, gptContext->lTrunc);
		gptContext->yTrunc[lGamma].gamma    = gamma;
		mkGptTblLine(& (gptContext->yTrunc[lGamma]) );
	}
}

/* Calculate starting and ending points for integral */
int calIntgLeng(int leng, double gamma, double lTrunc) {
	int intgLeng;
	if (lTrunc < EPS) { /* If lTrunc is equal zero, use full length */
		intgLeng = leng - 1;
	} else { /* Length of integral is decided according to gamma */
		intgLeng = (int) floor(lTrunc / sqrt(2.0 * gamma) + 0.5);
		intgLeng = (intgLeng < leng - 1) ? intgLeng : (leng - 1);
	}
#ifdef DEBUG
	printf("intgLeng = %d \n", intgLeng);
#endif
	return intgLeng;
}

/* Make tables of starting and ending points for integral */
void mkGptTblLine(GptTblLine *gptTblLine) {
	int    leng      = gptTblLine->leng;
	int    intgLeng  = gptTblLine->intgLeng;
	int    intgLeng2 = intgLeng * 2;
	double gamma     = gptTblLine->gamma;
	double sqGamma   = sqrt(gamma);
	int    st, gSt, endOver;
	int    pos;

	for (pos = -intgLeng ; pos <= intgLeng ; ++pos)
		gptTblLine->gauss[pos + intgLeng] = sqGamma * exp(- gamma * pos * pos);

	for (pos = 0 ; pos <= leng ; ++pos) {
		st  = pos - intgLeng;
		gSt = 0;
		if (st < 0) {
			gSt = -st;
			st = 0;
		}
		/* Basically pixels for integral is intgLeng * 2 + 1 */
		gptTblLine->st[pos]   = st;
		gptTblLine->gSt[pos]  = gSt;
		endOver               = (pos + intgLeng) - (leng - 1);
		gptTblLine->gEnd[pos] = (endOver > 0)? intgLeng2 - endOver : intgLeng2;
	}
}

/* Make templates for GPT calculation */
void mkGptTemplate(ImgInf *imgInf, GptContext *gptContext) {
	int nx            = gptContext->nx;
	int ny            = gptContext->ny;
	int nDirNd        = gptContext->nDirNd;
	int nxDirNd       = gptContext->nx * nDirNd;
	double *imgCan    = imgInf->can;
	int nxyDirNdGamma = gptContext->nxyDirNdGamma;
	int lGamma, offsetG;
	int x, y, pos, tPos, gPos, ynx, lDir;
	GptTblLine *xTblL, *yTblL;
	double *gh0, *gH0, *gh1, *gh2, *gH1x, *gH1y;
	double tx, ty, dx, dy, cx = gptContext->cx, cy = gptContext->cy;

	gh0   = gptContext->gh0;
	gh1   = gptContext->gh1;
	gh2   = gptContext->gh2;
	gH0   = imgInf->gH0;
	gH1x  = imgInf->gH1x;
	gH1y  = imgInf->gH1y;

	for (pos = 0 ; pos < nxyDirNdGamma ; ++pos) {
		gh0[pos]   = 0.0; gh1[pos]   = 0.0; gh2[pos]   = 0.0;
		gH0[pos]   = 0.0; gH1x[pos]  = 0.0; gH1y[pos]  = 0.0;
	}

	offsetG = 0; ynx = 0;
	for (lGamma = 0 ; lGamma < gptContext->nGamma ; ++lGamma) {
		if (gptContext->isTrunc == NOGAUSSTRUNC) {
			xTblL = & (gptContext->xFull[lGamma]) ;
			yTblL = & (gptContext->yFull[lGamma]) ;
		} else {
			xTblL = & (gptContext->xTrunc[lGamma]) ;
			yTblL = & (gptContext->yTrunc[lGamma]) ;
		}
		pos = 0;
		for (y = 0 ; y < ny ; ++y) {
			dy = y - cy;
			for (x = 0 ; x < nx ; ++x) {
#ifdef DEBUG
				printf("(%d, %d) st = %d gSt = %d gEnd = %d\n",  x, y, yTblL->st[y], yTblL->gSt[y], yTblL->gEnd[y]);
#endif
				gPos = yTblL->gSt[y];
				ty   = dy * imgCan[pos];
#ifdef DIRMATCH
				tPos = offsetG + (yTblL->st[y] * nx + x) * nDirNd + imgInf->ang[pos];
#else
				tPos = offsetG + (yTblL->st[y] * nx + x);
#endif
				ty   = dy * imgCan[pos];
				while(gPos <= yTblL->gEnd[y]) {
					gh0[tPos] += yTblL->gauss[gPos] * imgCan[pos];
					gh1[tPos] += yTblL->gauss[gPos] * ty;
					++gPos; tPos += nxDirNd;
				}
				++pos;
			}
		}

		pos = offsetG; ynx = 0;
		for (y = 0 ; y < ny ; ++y) {
			for (x = 0 ; x < nx ; ++x) {
				dx = x - cx;
				for (lDir = 0 ; lDir < nDirNd ; ++lDir) {
					tPos = offsetG + (ynx + xTblL->st[x]) * nDirNd + lDir;
#ifdef DEBUG
					printf("(%d, %d) dir = %d tPos = %d st = %d gSt = %d gEnd = %d\n", x, y, lDir, tPos, xTblL->st[x], xTblL->gSt[x], xTblL->gEnd[x]);
#endif
					gPos = xTblL->gSt[x];
					tx   = dx * gh0[pos];
					while(gPos <= xTblL->gEnd[x]) {
						gH0[tPos]  += xTblL->gauss[gPos] * gh0[pos];
						gH1x[tPos] += xTblL->gauss[gPos] * tx;
						gH1y[tPos] += xTblL->gauss[gPos] * gh1[pos];
						++gPos; tPos += nDirNd;
					}
					++pos; /* 方向も含まれる */
				}
			}
			ynx += nx;
		}
		offsetG += gptContext->nxyDirNd; /* For next gamma */
	}
}

void gptTransformPoint(double *gpt, double *inP, double *outP) {
	int i, j;
	double sum;
	for(i = 0 ; i < 3 ; ++i) {
		sum = 0.0;
		for(j = 0 ; j < 3 ; ++j) {
			sum += gpt[i + j * 3] * inP[j];
		}
		outP[i] = sum;
	}
}

void gptTransformGpt(double *gpt1, double *gpt2, double *outGpt) {
	int i, j, k;
	double sum;
	for(i = 0 ; i < 3 ; ++i) {
		for(j = 0 ; j < 3 ; ++j) {
			sum = 0.0;
			for(k = 0 ; k < 3 ; ++k) {
				sum += gpt1[i + k * 3] * gpt2[k + j * 3];
			}
			outGpt[i + j * 3] = sum;
		}
	}
}

/* Projection transformation of the image by bilinear interpolation */
void gptInvTransformImage(double *gpt, int *inImg, int *outImg, GptContext *gptContext){
	int x, y, xOrg, yOrg, posOrg, nx = gptContext->nx;
	double xOrgFl, xOrgFrac, yOrgFl, yOrgFrac;
	double inVect[3], outVect[3];

	/* Output image generation by bilinear interpolation */
	inVect[2] = 1.0;
	for (y = 0 ; y < gptContext->ny ; y++) {
		inVect[1] = y - gptContext->cy;
		for (x = 0 ; x < nx ; x++) {
			inVect[0] = x - gptContext->cx;
			gptTransformPoint(gpt, inVect, outVect);
			xOrgFl = outVect[0] / outVect[2] + gptContext->cx;
			yOrgFl = outVect[1] / outVect[2] + gptContext->cy;
			xOrg   = (int) floor(xOrgFl);
			yOrg   = (int) floor(yOrgFl);
			xOrgFrac = xOrgFl - xOrg;
			yOrgFrac = yOrgFl - yOrg;
			if (xOrg >= 0 && xOrg + 1 < nx && yOrg >= 0 && yOrg + 1 < gptContext->ny) {
				posOrg = yOrg * nx + xOrg;
				outImg[y * nx + x] = (int) (
						(1.0 - yOrgFrac) * ((1.0 - xOrgFrac) * inImg[posOrg] + xOrgFrac * inImg[posOrg + 1])
						+ yOrgFrac * ((1.0 - xOrgFrac) * inImg[posOrg + nx] + xOrgFrac * inImg[posOrg + nx + 1]));
			} else {
				outImg[y * nx + x] = WHITE;
			}
		}
	}
}

/* Projection transformation of the image by bilinear interpolation */
void gptInvTransformImageD(double *gpt, double *inImg, double *outImg, GptContext *gptContext){
	int x, y, xOrg, yOrg, posOrg, nx = gptContext->nx;
	double xOrgFl, xOrgFrac, yOrgFl, yOrgFrac;
	double inVect[3], outVect[3];

	/* Output image generation by bilinear interpolation */
	inVect[2] = 1.0;
	for (y = 0 ; y < gptContext->ny ; y++) {
		inVect[1] = y - gptContext->cy;
		for (x = 0 ; x < nx ; x++) {
			inVect[0] = x - gptContext->cx;
			gptTransformPoint(gpt, inVect, outVect);
			xOrgFl = outVect[0] / outVect[2] + gptContext->cx;
			yOrgFl = outVect[1] / outVect[2] + gptContext->cy;
			xOrg   = (int) floor(xOrgFl);
			yOrg   = (int) floor(yOrgFl);
			xOrgFrac = xOrgFl - xOrg;
			yOrgFrac = yOrgFl - yOrg;
			if (xOrg >= 0 && xOrg + 1 < nx && yOrg >= 0 && yOrg + 1 < gptContext->ny) {
				posOrg = yOrg * nx + xOrg;
				outImg[y * nx + x] = 	(1.0 - yOrgFrac) * ((1.0 - xOrgFrac) * inImg[posOrg] + xOrgFrac * inImg[posOrg + 1])
								+ yOrgFrac * ((1.0 - xOrgFrac) * inImg[posOrg + nx] + xOrgFrac * inImg[posOrg + nx + 1]);
			} else {
				outImg[y * nx + x] = 0.0;
			}
		}
	}
}

void gptInit(double *gpt) {
	gpt[C11] = 1.0; gpt[C12] = 0.0; gpt[C21] = 0.0; gpt[C22] = 1.0;
	gpt[C13] = 0.0; gpt[C23] = 0.0;
	gpt[C31] = 0.0; gpt[C32] = 0.0;
	gpt[C33] = 1.0;
}

/* Copy GPT */
void gptCopy(double *inGpt, double *outGpt) {
	int i;
	for(i = 0 ; i < 9 ; ++i) {
		outGpt[i] = inGpt[i];
	}
}

/* Inverse projection transformation */
int gptInverse(double *gpt, double *iGpt) {
	double det;
	det =   gpt[C11] * gpt[C22] * gpt[C33]
	                                  + gpt[C21] * gpt[C32] * gpt[C13]
	                                                              + gpt[C31] * gpt[C12] * gpt[C23]
	                                                                                          - gpt[C11] * gpt[C32] * gpt[C23]
	                                                                                                                      - gpt[C21] * gpt[C12] * gpt[C33]
	                                                                                                                                                  - gpt[C31] * gpt[C22] * gpt[C13];
	if (fabs(det) < EPS) {
		printf("Singular GPT is generated!!! in gptInverse() \n");
		gptPr(gpt, "Values");
		return 1;
	}
	iGpt[C11] =  (gpt[C22] * gpt[C33] - gpt[C32] * gpt[C23]) / det;
	iGpt[C21] = -(gpt[C21] * gpt[C33] - gpt[C31] * gpt[C23]) / det;
	iGpt[C31] =  (gpt[C21] * gpt[C32] - gpt[C31] * gpt[C22]) / det;

	iGpt[C12] = -(gpt[C12] * gpt[C33] - gpt[C32] * gpt[C13]) / det;
	iGpt[C22] =  (gpt[C11] * gpt[C33] - gpt[C31] * gpt[C13]) / det;
	iGpt[C32] = -(gpt[C11] * gpt[C32] - gpt[C31] * gpt[C12]) / det;

	iGpt[C13] =  (gpt[C12] * gpt[C23] - gpt[C22] * gpt[C13]) / det;
	iGpt[C23] = -(gpt[C11] * gpt[C23] - gpt[C21] * gpt[C13]) / det;
	iGpt[C33] =  (gpt[C11] * gpt[C22] - gpt[C21] * gpt[C12]) / det;
	return 0;
}


/* Display GPT calculation templates */
void gptTblPr(GptContext *context) {
	int x, y, pos, lGamma;
	GptTblLine *gptTblL;

	for (lGamma = 0 ; lGamma < context->nGamma ; ++lGamma) {
		gptTblL = & (context->yFull[lGamma]);
		for (y = 0 ; y < context->ny ; ++y) {
			printf("yFull dir = %d y = %d st = %5d gSt = %d gEnd = %d\n",
					lGamma, y, gptTblL->st[y], gptTblL->gSt[y], gptTblL->gEnd[y]);
		}
		for (pos = 0 ; pos <= 2 * gptTblL->intgLeng ; ++pos) {
			printf("pos = %d  gauss = %f \n", pos, gptTblL->gauss[pos]);
		}
	}

	for (lGamma = 0 ; lGamma < context->nGamma ; ++lGamma) {
		gptTblL = & (context->yTrunc[lGamma]);
		for (y = 0 ; y < context->ny ; ++y) {
			printf("yTrunc dir = %d y = %d st = %5d gSt = %d gEnd = %d\n",
					lGamma, y, gptTblL->st[y], gptTblL->gSt[y], gptTblL->gEnd[y]);
		}
		for (pos = 0 ; pos <= 2 * gptTblL->intgLeng ; ++pos) {
			printf("pos = %d  gauss = %f \n", pos, gptTblL->gauss[pos]);
		}
	}

	for (lGamma = 0 ; lGamma < context->nGamma ; ++lGamma) {
		gptTblL = & (context->xFull[lGamma]);
		for (x = 0 ; x < context->nx ; ++x) {
			printf("xFull dir = %d x = %d st = %5d gSt = %d gEnd = %d\n",
					lGamma, x, gptTblL->st[x], gptTblL->gSt[x], gptTblL->gEnd[x]);
		}
		for (pos = 0 ; pos <= 2 * gptTblL->intgLeng ; ++pos) {
			printf("pos = %d  gauss = %f \n", pos, gptTblL->gauss[pos]);
		}
	}

	for (lGamma = 0 ; lGamma < context->nGamma ; ++lGamma) {
		gptTblL = & (context->xTrunc[lGamma]);
		for (x = 0 ; x < context->nx ; ++x) {
			printf("xTrunc dir = %d x = %d st = %5d gSt = %d gEnd = %d\n",
					lGamma, x, gptTblL->st[x], gptTblL->gSt[x], gptTblL->gEnd[x]);
		}
		for (pos = 0 ; pos <= 2 * gptTblL->intgLeng ; ++pos) {
			printf("pos = %d  gauss = %f \n", pos, gptTblL->gauss[pos]);
		}
	}
}

/* Display GPT template */
void gptTemplatePr(GptContext *context) {
	int x, y, lDir, lGamma, pos;
	int nx = context->nx, ny = context->ny, nDirG = context->nDirNd;

	for (lGamma = 0 ; lGamma < context->nGamma ; ++lGamma) {
		for (lDir= 0 ; lDir < context->nDirNd ; ++lDir) {
			for (y = 0 ; y < context->ny ; ++y) {
				for (x = 0 ; x < context->nx ; ++x) {
					pos = lDir + nDirG * (x + nx * (y + ny * lGamma));
					printf("lGamma = %d (%3d, %3d) dir =%d %12.7f  %12.7f   %12.7f  %12.7f  %12.7f\n",
							lGamma, x, y, lDir, context->gh0[pos], context->gh1[pos],
							context->tgtImgInf->gH0[pos], context->tgtImgInf->gH1x[pos], context->tgtImgInf->gH1y[pos]);
					++pos;
				}
			}
		}
	}
}
