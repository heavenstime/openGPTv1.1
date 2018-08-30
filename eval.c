/*
 * eval.c
 *
 *  Created on: 2013/07/02
 *      Author: yamasita
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "eval.h"
#include "gpt.h"
#include "parameters.h"
#include "utility.h"

//#define DEBUG
//#define L2  /* correlation or L2 difference */
//#define DEBUGWDCH

double nndegdP0(int xOrg, int yOrg, int angCode, int *imgAng, NndegdContext *nndegdContext);
double nndegdP1(int xOrg, int yOrg, int posOrg, int angCode, int *imgAng, NndegdContext *nndegdContext);
void wdchCal(int *ang, double *hist4, int nx, int ny);

/* Calculate image correlation */
double imgCor(double *img1, double *img2, int nxy) {
	int pos;
	double cor = 0.0;
#ifdef L2
	for (pos = 0 ; pos < nxy ; ++pos) cor -= (img1[pos] - img2[pos]) * (img1[pos] - img2[pos]);
#else
	for (pos = 0 ; pos < nxy ; ++pos) cor += img1[pos] * img2[pos];
#endif
	return cor;
}

/* Allocate region for calculation of Dnn */
void nndegdInit(int nx, int ny, NndegdContext *nndegdContext) {
	/* Parameter independent of content of image */
	nndegdContext->nx          = nx = COL;
	nndegdContext->ny          = ny = ROW;
	nndegdContext->nDir        = DIRECTION;
	/* Constant decided from parameters */
	nndegdContext->nst           = (2 * nx - 1) * (2 * ny - 1);
	nndegdContext->dDiagSq       = nx * nx + ny * ny;
	nndegdContext->dDiag         = sqrt((double) nndegdContext->dDiagSq);
	/* Tables which do not depend on content of image */
	nndegdContext->nndegdBound  = (int *) malloc(sizeof(int) * nx * ny);
	nndegdContext->scanNodes    = (ScanNode *) malloc(sizeof(ScanNode) * nndegdContext->nst);
}

/* Calculate weight according to curvature */
void calWeight(int *img, double *weight, int nx, int ny, double *workImgD1, double *workImgD2, double *workImgD3, double *workImgD4) {
	int nIter = CURVITER;
	int ix, iy, pos, posGrad;  /* Pos variable */
	int nx1  = nx + 1;
	int ny1  = ny + 1;
	int nxy1 = nx1 * ny1;
	int curvSmooth = - CURVSMOOTH;
	double delta_RD, delta_LD, tr, det;

	int  pos33[9] = {-nx1 - 1, -nx1, -nx1 + 1, -1, 0, 1, nx1 + 1, nx1, nx1 - 1};
	double w33[9] = BLUERCOEFLIST; /* Approximated Gaussian weight */

	imgInfSmoothD(img, workImgD4, curvSmooth, workImgD1, nx, ny);
	/* Calulate differentials */
	for (posGrad = 0 ; posGrad < nxy1 ; ++posGrad)
		workImgD1[posGrad] = workImgD2[posGrad] = workImgD3[posGrad] = 0.0;

	pos     = 0;
	posGrad = nx1 + 1;
	for (iy = 0 ; iy < ny - 1 ; ++iy) {
		for (ix = 0 ; ix < nx - 1 ; ++ix) {
			delta_RD = workImgD4[pos + 1] - workImgD4[pos + nx];
			delta_LD = workImgD4[pos]     - workImgD4[pos + nx + 1];
			workImgD1[posGrad] = delta_RD * delta_RD;
			workImgD2[posGrad] = delta_RD * delta_LD;
			workImgD3[posGrad] = delta_LD * delta_LD;
			++pos; ++posGrad;
		}
		++pos; /* x < nx - 1 のため */
		posGrad += 2;
	}
#ifdef DEBUGWEIGHT
	printf("Image data\n");
	pos = 0;
	for (iy = 0 ; iy < ny  ; ++iy) {
		for (ix = 0 ; ix < nx ; ++ix) {
			printf("%5d ", img[pos++]);
		}
		printf("\n");
	}
	printf("\nGrad1 Grad1\n");

	pos = 0;
	for (iy = 0 ; iy < ny1  ; ++iy) {
		for (ix = 0 ; ix < nx1 ; ++ix) {
			printf("%5.1f ", workImgD1[pos++]);
		}
		printf("\n");
	}
	printf("\n Grad1 Grad2\n");
	pos = 0;
	for (iy = 0 ; iy < ny1  ; ++iy) {
		for (ix = 0 ; ix < nx1 ; ++ix) {
			printf("%5.1f ", workImgD2[pos++]);
		}
		printf("\n");
	}
	printf("\n Grad2 Grad2\n");
	pos = 0;
	for (iy = 0 ; iy < ny1  ; ++iy) {
		for (ix = 0 ; ix < nx1 ; ++ix) {
			printf("%5.1f ", workImgD3[pos++]);
		}
		printf("\n");
	}
	printf("\n");
#endif

	for (posGrad = 0 ; posGrad < nxy1 ; ++posGrad) workImgD4[posGrad] = 0.0;
	bluer33Iter(workImgD1, workImgD4, pos33, w33, nIter, nx - 1, ny - 1);
	bluer33Iter(workImgD2, workImgD4, pos33, w33, nIter, nx - 1, ny - 1);
	bluer33Iter(workImgD3, workImgD4, pos33, w33, nIter, nx - 1, ny - 1);

	pos = 0; posGrad = nx1 + 1;
	for (iy = 0 ; iy < ny - 1 ; ++iy) {
		for (ix = 0 ; ix < nx - 1 ; ++ix) {
			tr  = workImgD1[posGrad] + workImgD3[posGrad];
			det = workImgD1[posGrad] * workImgD3[posGrad] - workImgD2[posGrad] * workImgD2[posGrad];
			// printf("(%d %d) = %f  %f %f\n", ix, iy, tr, det, det / (tr * tr));
			weight[pos] = CURVCOEFF * ((tr > EPS)? det / (tr * tr) : 0.0) + 1.0;
			++pos; ++posGrad;
		}
		weight[pos++] = 1.0; /* x < nx - 1 のため */
		posGrad += 2;
	}
	for (ix = 0 ; ix < nx ; ++ix) weight[pos++] = 1.0; /* for iy = ny - 1 */
}

/* calculation of mean of nearest-neighbor interpoint distances */
/* with the same angle code between two images */
double nndegd(int nndegdType, int *imgAng1, int *imgAng2, NndegdContext *nndegdContext, double *nndegdTbl) {
	/* Average of both side */
	int x, y, pos, count1 = 0, count2 = 0, tPos;
	double nndegd1 = 0.0, nndegd2 = 0.0;
	int nDir = nndegdContext->nDir;

	switch (nndegdType) {
	case NNDEGDORG:
		pos = 0;
		for (y = 0 ; y < nndegdContext->ny ; ++y) {
			for (x = 0 ; x < nndegdContext->nx ; ++x) {
				if (imgAng1[pos] != NoDIRECTION) {
					nndegd1 += nndegdP0(x, y, imgAng1[pos], imgAng2, nndegdContext); ++count1;
				}
				if (imgAng2[pos] != NoDIRECTION) {
					nndegd2 += nndegdP0(x, y, imgAng2[pos], imgAng1, nndegdContext); ++count2;
				}
				++pos;
			}
		}
		break;
	case NNDEGDTBL:
		pos = 0;
		for (y = 0 ; y < nndegdContext->ny ; ++y) {
			for (x = 0 ; x < nndegdContext->nx ; ++x) {
				if (imgAng1[pos] != NoDIRECTION) {
					nndegd1 += nndegdP1(x, y, pos, imgAng1[pos], imgAng2, nndegdContext); ++count1;
				}
				if (imgAng2[pos] != NoDIRECTION) {
					nndegd2 += nndegdP1(x, y, pos, imgAng2[pos], imgAng1, nndegdContext); ++count2;
				}
				++pos;
			}
		}
		break;
	case NNDEGDTMPLT:
		pos = 0; tPos = 0;
		for (y = 0 ; y < nndegdContext->ny ; ++y) {
			for (x = 0 ; x < nndegdContext->nx ; ++x) {
				if (imgAng1[pos] != NoDIRECTION) {
					nndegd1 += nndegdTbl[tPos + imgAng1[pos]]; ++count1;
				}
				if (imgAng2[pos] != NoDIRECTION) {
					nndegd2 += nndegdP1(x, y, pos, imgAng2[pos], imgAng1, nndegdContext); ++count2;
				}
				++pos; tPos += nDir;
			}
		}
		break;
	case NNDEGDTMPLTHALF: /* Omit calculation of nndegd2 (= nndegd1) */
		pos = 0; tPos = 0;
		for (y = 0 ; y < nndegdContext->ny ; ++y) {
			for (x = 0 ; x < nndegdContext->nx ; ++x) {
				if (imgAng1[pos] != NoDIRECTION) {
					nndegd1 += nndegdTbl[tPos + imgAng1[pos]]; ++count1;
				}
				++pos; tPos += nDir;
			}
		}
		nndegd2   =  nndegd1; /* nndegd2 is set according to nndegd1 */
		count2 = count1;
		break;
	}

#ifdef DEBUG
	if (count1 != 0) 	printf("nndegd  nndegd1 = %f ", nndegd1 / (double) count1);
	else 	printf("nndegd  count1 = 0 ");
	if (count2 != 0) 	printf(" nndegd2 = %f\n", nndegd2 / (double) count2);
	else 	printf("count2 = 0\n");
#endif
	if (count1 <= 1 || count2 <= 1) {
		if (count1 > 1) 	printf("nndegd  nndegd1 = %f ", nndegd1 / (double) count1);
		else 	printf("nndegd  count1 <= 1 ");
		if (count2 > 1) 	printf(" nndegd2 = %f\n", nndegd2 / (double) count2);
		else 	printf("count2 <= 1\n");
		return NNDEGDFAIL;
	}
	return 0.5 * (nndegd1 / (double) count1 + nndegd2 / (double) count2);
}

double nndegdP0(int xOrg, int yOrg, int angCode, int *imgAng, NndegdContext *nndegdContext) {
	int  x, y, nx = nndegdContext->nx, pos = 0;
	int  delta, min;
	min = nndegdContext->dDiagSq;
	for (y = 0 ; y < nndegdContext->ny ; ++y) {
		for (x = 0 ; x < nx ; ++x) {
			if (imgAng[pos] == angCode) {
				delta = (y - yOrg) * (y - yOrg) + (x - xOrg) * (x - xOrg);
				if (delta < min) min = delta;
			}
			++pos;
		}
	}
	return sqrt((double) min);
}

/* Calculation of mean of nearest-neighbor interpoint distances */
/* with the same angle code between two images using ScanTbl */
double nndegdP1(int xOrg, int yOrg, int posOrg, int angCode, int *imgAng, NndegdContext *nndegdContext) {
	int  sx, sy, sPos, dBound;
	ScanNode *sNodes;

	/* Initialize */
	sNodes = nndegdContext->scanNodes;
	dBound = (nndegdContext->nndegdBound)[posOrg];
	for (sPos = 0 ; sPos < nndegdContext->nst; ++sPos) {
		if (sNodes[sPos].l1Dist > dBound) break;
		if (imgAng[posOrg + sNodes[sPos].dPos] == angCode) {
			return sNodes[sPos].dist;
		}
	}
	for ( ; sPos < nndegdContext->nst; ++sPos) {
		sy = yOrg + sNodes[sPos].dy;
		sx = xOrg + sNodes[sPos].dx;
		if (sx >= 0 && sy >= 0 && sx < nndegdContext->nx && sy < nndegdContext->ny) {
			if (imgAng[posOrg + sNodes[sPos].dPos] == angCode) {
				return sNodes[sPos].dist;
			}
		}
	}
	return (nndegdContext->dDiag);
}

/* Make scan table */
void mkScanTbl(NndegdContext *nndegdContext) {
	int x, y, nx = nndegdContext->nx, ny = nndegdContext->ny, nPos, sPos, tPos, max, min;
	ScanNode *sNodes, *minNode, tmpNode;
	int *nndegdBound;

	/* Initialize */
	sNodes      = nndegdContext->scanNodes;
	nndegdBound = nndegdContext->nndegdBound;

	sPos = 0;
	for (y = - ny + 1 ; y < ny ; ++y) {
		for (x = - nx + 1 ; x < nx ; ++x) {
			sNodes[sPos].dx = x;
			sNodes[sPos].dy = y;
			sNodes[sPos].dist = sqrt((double) (x * x + y * y));
			sNodes[sPos].dPos = y * nx + x;
			max = abs(x);
			if (max < abs(y)) max = abs(y);
			sNodes[sPos].l1Dist = max;
#ifdef DEBUGSCANTBL
			printf("(%d, %d) sPos = %d liDist = %d dist = %f \n", sNodes[sPos].dx, sNodes[sPos].dy, sNodes[sPos].dPos, sNodes[sPos].l1Dist, sNodes[sPos].dist);
#endif
			++sPos;
		}
	}
	nPos = sPos;

	/* Sort by distance in ascending order */
	for (sPos = 0 ; sPos < nPos ; ++sPos) {
		minNode = & sNodes[sPos];
		for (tPos = sPos ; tPos < nPos ; ++tPos) {
			if (minNode->dist > sNodes[tPos].dist) minNode = & sNodes[tPos];
		}
		tmpNode      = sNodes[sPos];
		sNodes[sPos] = *minNode;
		*minNode     = tmpNode;
	}

	/* make bound of distance (including square) */
	for (y = 0 ; y < ny ; ++y) {
		for (x = 0 ; x < nx ; ++x) {
			min = x;
			if (min > y) min = y;
			if (min > nx - x - 1) min = nx - x - 1;
			if (min > ny - y - 1) min = ny - y - 1;
			nndegdBound[y * nx + x] = min;
		}
	}
}

/* Make table of distance of the same direction */
void mkNnDistTbl(int *imgAng, double *nndegdTbl, NndegdContext *nndegdContext) {
	int x, y, pos, nDir = nndegdContext->nDir;
	int tblPos, finishAll, finish[nDir], i;
	int sx, sy, sPos, angCode, dBound;
	ScanNode *sNodes;

	if (nndegdContext->nndegdType < 2) return;

	/* Initialize */
	sNodes    = nndegdContext->scanNodes;
	for (pos = 0 ; pos < nndegdContext->nx * nndegdContext->ny * nDir; ++pos) {
		nndegdTbl[pos] = nndegdContext->dDiag;
	}
	pos = 0; tblPos = 0;
	for (y = 0 ; y < nndegdContext->ny ; ++y) {
		for (x = 0 ; x < nndegdContext->nx ; ++x) {
			dBound = nndegdContext->nndegdBound[pos];
			finishAll = nDir;
			for (i = 0 ; i < nDir; ++i) finish[i] = 0;
			for (sPos = 0 ; sPos < nndegdContext->nst; ++sPos) {
#ifdef DEBUG
				printf("A pos = %d sPos = %d \n", pos, sPos);
#endif
				if (sNodes[sPos].l1Dist > dBound || finishAll == 0) break;
				angCode = imgAng[pos + sNodes[sPos].dPos];
				if(angCode != NoDIRECTION && finish[angCode] == 0) {
#ifdef DEBUG
					printf("ang A %d  %d\n", angCode, pos + sNodes[sPos].dPos);
#endif
					nndegdTbl[tblPos + angCode] = sNodes[sPos].dist;
					finish[angCode] = 1;
					--finishAll;
				}
			}
			for ( ; sPos < nndegdContext->nst; ++sPos) {
				if (finishAll == 0) break; /* Distances are obtained with all directions */
				sy = y + sNodes[sPos].dy;
				sx = x + sNodes[sPos].dx;
				if (sx >= 0 && sy >= 0 && sx < nndegdContext->nx && sy < nndegdContext->ny) {
					angCode = imgAng[pos + sNodes[sPos].dPos];
					if(angCode != NoDIRECTION && finish[angCode] == 0) {
#ifdef DEBUG
						printf("ang B %d  %d\n", angCode, pos + sNodes[sPos].dPos);
#endif
						nndegdTbl[tblPos + angCode] = sNodes[sPos].dist;
						finish[angCode] = 1;
						--finishAll;
					}
				}
			}
			++pos; tblPos += nDir;
		}
	}
}

/* evType = 1 - 4 */
double nndegdEv(int nndegdType, int evType, int *imgAng1, int *imgAng2, double *weight1, double *weight2, NndegdContext *nndegdContext, double *nndegdTbl) {
	/* Average of both side */
	int ix, iy, pos, tPos;
	double w1, w2, totalW1 = 0.0, totalW2 = 0.0;
	double nndegdPnt1, nndegdPnt2, nndegd1 = 0.0, nndegd2 = 0.0, nndegdSq1 = 0.0, nndegdSq2 = 0.0, val = 0.0;
	int nDir = nndegdContext->nDir;

	pos = 0; tPos = 0;
	for (iy = 0 ; iy < nndegdContext->ny ; ++iy) {
		for (ix = 0 ; ix < nndegdContext->nx ; ++ix) {
			nndegdPnt1 = nndegdPnt2 = 0.0; w1 = w2 = 0.0;
			/* For forward measure */
			if (imgAng1[pos] != NoDIRECTION) {
				nndegdPnt1 = nndegdTbl[tPos + imgAng1[pos]];
				switch(evType) {
				case EVALNNDEGD:
				case EVALNNDEGDVAR:
				case EVALNNDEGDSQR:
				case EVALNNDEGDSQRT:
					w1 = 1.0; break;
				case EVALNNDEGDWEIGHT:
				case EVALNNDEGDSQRWEIGHT:
					w1 = weight1[pos]; break;
				}
			}

			/* For backward measure */
			switch (nndegdType) {
			case NNDEGDTMPLT:
				if (imgAng2[pos] != NoDIRECTION) {
					nndegdPnt2 = nndegdP1(ix, iy, pos, imgAng2[pos], imgAng1, nndegdContext);
					switch(evType) {
					case EVALNNDEGD:
					case EVALNNDEGDVAR:
					case EVALNNDEGDSQR:
					case EVALNNDEGDSQRT:
						w2 = 1.0; break;
					case EVALNNDEGDWEIGHT:
					case EVALNNDEGDSQRWEIGHT:
						w2 = weight2[pos]; break;
					}
				}
				break;
			case NNDEGDTMPLTHALF: /* Omit calculation of nndegd2 (= nndegd1) */
				nndegdPnt2 = nndegdPnt1;
				w2    = w1;
				break;
			}

			switch(evType) {
			case EVALNNDEGD:
				nndegd1 += nndegdPnt1; nndegd2 += nndegdPnt2; break;
			case EVALNNDEGDVAR:
				nndegd1 += nndegdPnt1; nndegdSq1 += nndegdPnt1 * nndegdPnt1;
				nndegd2 += nndegdPnt2; nndegdSq2 += nndegdPnt2 * nndegdPnt2; break;
			case EVALNNDEGDSQR:
				nndegdSq1 += nndegdPnt1 * nndegdPnt1; nndegdSq2 += nndegdPnt2 * nndegdPnt2; break;
			case EVALNNDEGDSQRT:
				nndegdSq1 += sqrt(nndegdPnt1); nndegdSq2 += sqrt(nndegdPnt2); break;
			case EVALNNDEGDWEIGHT:
				nndegd1 += w1 * nndegdPnt1; nndegd2 += w2 * nndegdPnt2; break;
			case EVALNNDEGDSQRWEIGHT:
				nndegdSq1 += w1 * nndegdPnt1 * nndegdPnt1; nndegdSq2 += w2 * nndegdPnt2 * nndegdPnt2; break;
			}
			totalW1 += w1; totalW2 += w2;
			//printf("ang1 = %d ang2 = %d nndegd1 = %f nndegd2 = %f count1 = %d count2 = %d \n", imgAng2[pos], imgAng2[pos], nndegd1, nndegd2, count1, count2);
			++pos; tPos += nDir;
		}
	}

	if (totalW1 <= 1.01 || totalW2 <= 1.01) {
		if (totalW1 > 1) 	printf("nndegdSq  nndegd1Sq = %f ", nndegd1 / (double) totalW1);
		else 	printf("nndegdSq  count1 <= 1 ");
		if (totalW2 > 1) 	printf(" nndegdSq2 = %f\n", nndegd2 / (double) totalW2);
		else 	printf("count2Sq <= 1\n");
		return NNDEGDFAIL;
	}
	switch(evType) {
	case EVALNNDEGD:
	case EVALNNDEGDWEIGHT:
		val = 0.5 * (nndegd1 / (double) totalW1 + nndegd2 / (double) totalW2); break;
	case EVALNNDEGDVAR:
		val = 0.5 * ((nndegdSq1 - nndegd1 * nndegd1 / (double) totalW1) / (double) (totalW1 - 1)
				+ (nndegdSq2 - nndegd2 * nndegd2 / (double) totalW2) / (double) (totalW2 - 1) );
		break;
	case EVALNNDEGDSQR:
	case EVALNNDEGDSQRT:
	case EVALNNDEGDSQRWEIGHT:
		val = 0.5 * (nndegdSq1 / (double) totalW1 + nndegdSq2 / (double) totalW2); break;
	}
	return val;
}

/* Calculate image difference with erosion */
double diffImgEr(int *img1, int *img2, int nErosion, int nx, int ny, double *workImgD1, double *workImgD2, double *workImgD3, double *workImgD4, int *diffErImg) {
	double coef = 0.1, diff, diff1 = 0.0, diff2 = 0.0;
	double mean, var, var1L1, var2L1;
	int dilerType = 5; /* Now 0 is background */
	int ixy, nxy = nx * ny;;

	for (ixy = 0 ; ixy <	nxy ; ++ixy) {
		diff = img1[ixy] - img2[ixy];
		if (diff >= 0) {
			workImgD1[ixy] = diff;
			workImgD4[ixy] = 0.0;
		} else {
			workImgD1[ixy] = 0.0;
			workImgD4[ixy] = -diff;
		}
	}

	dilationIter(dilerType, workImgD1, workImgD2, nx, ny, coef, nErosion);
	for (ixy = 0 ; ixy < nxy ; ++ixy) {
		diff1          += workImgD2[ixy];
#ifdef DEBUGDIFFER
		diffErImg[ixy]  = workImgD2[ixy];
#endif
	}
	dilationIter(dilerType, workImgD4, workImgD2, nx, ny, coef, nErosion);
	for (ixy = 0 ; ixy < nxy ; ++ixy) {
		diff2          += workImgD2[ixy];
#ifdef DEBUGDIFFER
		diffErImg[ixy] += workImgD2[ixy];
#endif
	}
	diff = diff1 + diff2;
	// diff = (diff1 > diff2) ? diff1 : diff2;
	calMeanVar(img1, &mean, &var, &var1L1, nxy);
	calMeanVar(img2, &mean, &var, &var2L1, nxy);
	return diff / (var1L1 + var2L1);
}

/* Weighted directional histgram */
double wdchEv(int *ang1, int *ang2, int nx, int ny) {
	static double hist41[128];
	static double hist42[128];
	double eval = 0.0;
	int ip;
	int n4 = 4 * 4 * 8;

	wdchCal(ang1, hist41, nx, ny);
	wdchCal(ang2, hist42, nx, ny);

	for (ip = 0 ; ip < n4 ; ++ip) eval += hist41[ip] * hist42[ip];
return eval;
}

void wdchCal(int *ang, double *hist4, int nx, int ny) {
	static int hist7e[968];
	static double hist4t[128];
	static double wdchGauss[] = WDCHGAUSS;
	static double wdchDirSmooth[] = WDCHDIRSMOOTH;
	double sum, total;
	int ix, iy, iDir, iDirD, ip, hx, hy, pos;
	int nx7e = 11, ny7e = 11, nDir = 8;
	int n7e = nx7e * ny7e * nDir;
	int n4 = 4 * 4 * 8;

	for (ip = 0 ; ip < n7e ; ++ip) hist7e[ip] = 0.0;

	for (iy = 0 ; iy < ny ; ++iy)  {
		hy = (int) (7.0 * iy / ny) + 2;
		for (ix = 0 ; ix < nx ; ++ix)  {
			hx = (int) (7.0 * ix / nx) + 2;
			if (ang[ix + nx * iy] != NoDIRECTION) {
				pos = ang[ix + nx * iy] + nDir * (hx + nx7e * hy);
				hist7e[pos] += 1;
			}
		}
	}

	/* Spatial smoothing and subsampling */
	for (iy = 0 ; iy < 4 ; ++iy)  {
		for (ix = 0 ; ix < 4 ; ++ix)  {
			for (iDir = 0 ; iDir < nDir ; ++iDir)  {
				sum = 0.0;
				for (hy = 0 ; hy < 5 ; ++hy)  {
					for (hx = 0 ; hx < 5 ; ++hx)  {
						sum += wdchGauss[hx + 5 * hy] * hist7e[iDir + nDir * (hx + 2 * ix + nx7e * (hy + 2 * iy))];
					}
				}
				hist4t[iDir + nDir * (ix + 4 * iy)] = sum;
			}
		}
	}

	/* Directional smoothing */
	total = 0.0;
	for (iy = 0 ; iy < 4 ; ++iy)  {
		for (ix = 0 ; ix < 4 ; ++ix)  {
			for (iDir = 0 ; iDir < nDir ; ++iDir)  {
				sum = 0.0;
				for (iDirD = -1 ; iDirD <= 1 ; ++iDirD)  {
						sum += wdchDirSmooth[iDirD + 1] * hist4t[(iDir + iDirD + nDir) % nDir + nDir * (ix + 4 * iy)];
				}
				hist4[iDir + nDir * (ix + 4 * iy)] = sum;
				total += sum * sum;
			}
		}
	}

	/* Norm normaization */
	for (ip = 0 ; ip < n4 ; ++ip) hist4[ip] /= sqrt(total);

#ifdef DEBUGWDCH
	static int flag = 0;
	if (flag == 0) {
		flag = 1;
		pos = 0;
		printf("Image data\n");
		for (iy = 0 ; iy < ny  ; ++iy) {
			if (iy % 4 == 0) printf("\n");
			for (ix = 0 ; ix < nx ; ++ix) {
				if (ix % 4 == 0) printf("  ");
				printf("%5d ", ang[pos++]);
			}
			printf("\n");
		}
		for (iy = 0 ; iy < 7  ; ++iy) {
			for (ix = 0 ; ix < 7 ; ++ix) {
				printf("(%3d %3d)   ", ix, iy);
				for (iDir = 0 ; iDir < 8 ; ++iDir) {
					printf("%3d ", hist7e[iDir + nDir * (ix + 2 + nx7e * (iy + 2))]);
				}
				printf("\n");
			}
		}
		for (iy = 0 ; iy < 4  ; ++iy) {
			for (ix = 0 ; ix < 4 ; ++ix) {
				printf("(%3d %3d) hist4t  ", ix, iy);
				for (iDir = 0 ; iDir < 8 ; ++iDir) {
					printf("%3f ", hist4t[iDir + nDir * (ix + 4 * iy)]);
				}
				printf("\n");
			}
		}
		for (iy = 0 ; iy < 4  ; ++iy) {
			for (ix = 0 ; ix < 4 ; ++ix) {
				printf("(%3d %3d)  hist4  ", ix, iy);
				for (iDir = 0 ; iDir < 8 ; ++iDir) {
					printf("%3f ", hist4[iDir + nDir * (ix + 4 * iy)]);
				}
				printf("\n");
			}
		}
	}
#endif

}

#ifdef DEBUGBLUER
/* Debug bluer33 */
int main(int argc, char *argv[]) {
	int nx = 10,      ny = 15;
	int nx2 = nx + 2, ny2 = ny + 2;
	double workImgD1[nx2 * ny2], workImgD4[nx2 * ny2];
	int ix, iy;
	int  pos33[9] = {-nx2 - 1, -nx2, -nx2 + 1, -1, 0, 1, nx2 + 1, nx2, nx2 - 1};
	double w33[9] = BLUERCOEFLIST;
	double sum = 0.0;

	for(ix = 0 ; ix < nx2 * ny2 ; ++ix) {
		workImgD1[ix] = workImgD4[ix] = 0.0;
	}
	workImgD1[4 + 1 + nx2 * (4 + 1)] = 1.0;
	bluer33Iter(workImgD1, workImgD4, pos33, w33, 3, 10, 10);
	for (iy = 0 ; iy < ny2 ; ++iy) {
		for (ix = 0 ; ix < nx2 ; ++ix) {
			sum += workImgD1[ix + nx2 * iy];
			printf("%6.3f ", workImgD1[ix + nx2 * iy]);
		}
		printf("\n");
	}
	printf("Sum = %f", sum);
	return 0;
}
#endif
