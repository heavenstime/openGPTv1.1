/*
 * utility.c
 *
 *  Created on: 2012/11/02
 *      Author: Yamasihta and Wakahara
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eval.h"
#include "gpt.h"
#include "parameters.h"
#include "utility.h"

//#define DEBUGDILATE

void swap_int(int *array, int i, int j);
void dilation(int type, double *inImgD, double *outImgD, int nx, int ny, double coef);
double hingeFunc(double in);
void bluer33(double *extInImg, double *extOutImg, int *pos33, double *w33, int nx, int ny);

/* Input of header & body information of pgm file (nExt : extension width */
void load_image_file(char *fileName, int *img, int nx, int ny, int nExt){
	/* unsigned char buffer[MAX_BUFFERSIZE];*/
	char buffer[MAX_BUFFERSIZE];
	FILE *fp;         /* File pointer */
	int ntx, nty;     /* Original image size */
	int max_gray;    /* Maximum gray level */
	int pos, ix, iy; /* Loop variable */

	/* Input file open */
	//printf("Input file = %s\n", fileName);
	fp = fopen(fileName, "rb");
	if (NULL == fp) {
		printf("     The file doesn't exist! %s \n\n", fileName);
		exit(1);
	}

	/* Check of file-type ---P5 */
	if (fgets(buffer, MAX_BUFFERSIZE, fp) == 0) {
		printf("     Read error of image file %s \n\n", fileName);
		exit(1);
	}
	if (buffer[0] != 'P' || buffer[1] != '5') {
		printf("     Mistaken file format, not P5! %s \n\n", fileName);
		exit(1);
	}

	/* input of x_size1, y_size1 */
	ntx = 0; nty = 0;
	while (ntx == 0 || nty == 0) {
		if (fgets(buffer, MAX_BUFFERSIZE, fp) == 0) {
			printf("     Read error of image file %s \n\n", fileName);
			exit(1);
		}
		if (buffer[0] != '#')  sscanf(buffer, "%d %d", & ntx, & nty);
	}

	if (ntx != (nx - 2 * nExt) || nty != (ny - 2 * nExt)) {
		printf("     Image size is not correct ! %s\n\n", fileName);
		exit(1);
	}

	/* input of max_gray */
	max_gray = 0;
	while (max_gray == 0) {
		if (fgets(buffer, MAX_BUFFERSIZE, fp) == 0) {
			printf("     Read error of image file %s \n\n", fileName);
			exit(1);
		}
		if (buffer[0] != '#')  sscanf(buffer, "%d", &max_gray);
	}

	if (max_gray != WHITE) {
		printf("     Invalid value of maximum gray level! %s \n\n", fileName);
		exit(1);
	}

	/* Clear buffer */
	for (pos = 0 ; pos < nx * ny ; ++pos) img[pos] = WHITE;
	/* Input of image data */
	for (iy = 0 ; iy < nty ; ++iy) {
		pos = nExt + (iy + nExt) * nx;
		for (ix = 0 ; ix < ntx ; ++ix) {
			img[pos++] = (int) ((unsigned char) fgetc(fp));
		}
	}
	fclose(fp);
}

/* Output of *img, nx, ny */
void save_image_file(char *filename, int *img, int nx, int ny){
	FILE *fp; /* File pointer */
	int  pos; /* Loop variable */

	fp = fopen(filename, "wb");
	/* output of pgm file header information */
	fputs("P5\n", fp);
	fputs("# Created by Image Processing\n", fp);
	fprintf(fp, "%d %d\n", nx, ny);
	fprintf(fp, "%d\n", WHITE);

	/* Output of image data */
	for (pos = 0 ; pos < nx * ny ; ++pos) {
		if(img[pos] > 255) fputc(WHITE, fp);
		else if(img[pos] < 0) fputc(BLACK, fp);
		else fputc(img[pos], fp);
	}
	fclose(fp);
}

void allcImgInf(ImgInf *imgInf, GptContext *gptContext) {
	imgInf->img       = (int *) malloc(sizeof(int) * gptContext->nxy);
	imgInf->ang       = (int *) malloc(sizeof(int) * gptContext->nxy);
	imgInf->can       = (double *) malloc(sizeof(double) * gptContext->nxy);
	imgInf->weight    = (double *) malloc(sizeof(double) * gptContext->nxy);
#ifndef NOGPT
	imgInf->gH0       = (double *) malloc(sizeof(double) * gptContext->nxyDirNdGamma);
	imgInf->gH1x      = (double *) malloc(sizeof(double) * gptContext->nxyDirNdGamma);
	imgInf->gH1y      = (double *) malloc(sizeof(double) * gptContext->nxyDirNdGamma);
#endif
	imgInf->nndegdTbl = (double *) malloc(sizeof(double) * gptContext->nxyDir);
	/* For dilation/erosion */
	if (gptContext->nDiler != 1) {
		imgInf->imgDiler       = (int *) malloc(sizeof(int) * gptContext->nxy * (gptContext->nDiler - 1));
		imgInf->angDiler       = (int *) malloc(sizeof(int) * gptContext->nxy * (gptContext->nDiler - 1));
		imgInf->canDiler       = (double *) malloc(sizeof(double) * gptContext->nxy * (gptContext->nDiler - 1));
		imgInf->weightDiler    = (double *) malloc(sizeof(double) * gptContext->nxy * (gptContext->nDiler - 1));
		imgInf->nndegdTblDiler = (double *) malloc(sizeof(double) * gptContext->nxyDir  * (gptContext->nDiler - 1));
	}
}

/* Read image and make data for GPT */
void loadMkImgInf(char *fileName, ImgInf *imgInf, GptContext *gptContext) {
	int dilerN, ixy, offset, nxy = gptContext->nxy;
	/* Read image */
	load_image_file(fileName, imgInf->img, gptContext->nx, gptContext->ny, gptContext->nExt);
	/* Make input image smooth by diferential morphology */
	if (gptContext->smooth != 0) imgInfSmooth(imgInf->img, - gptContext->smooth, gptContext);

	/* Canonical Definition */
	defcan(imgInf->img, imgInf->can, gptContext->nxy);
	/* Make an edge direction image */
	roberts8(imgInf->img, imgInf->ang, gptContext->nx, gptContext->ny);
	/* Make template for distance of the same direction */
	mkNnDistTbl(imgInf->ang, imgInf->nndegdTbl, gptContext->nndegdContext);
	/* Make GPT template */
#ifndef NOGPT
	mkGptTemplate(imgInf, gptContext);
#endif
	/* Make weight data */
	if (gptContext->evalType == EVALNNDEGDWEIGHT || gptContext->evalType == EVALNNDEGDSQRWEIGHT)
		calWeight(imgInf->img, imgInf->weight, gptContext->nx, gptContext->ny, gptContext->workImgD1, gptContext->workImgD2, gptContext->workImgD3, gptContext->workImgD4);

	/* for dilation/Erosion data */
	for (dilerN = 1 ; dilerN < gptContext->nDiler ; ++dilerN) {
		offset = (dilerN - 1) * nxy;
		for (ixy = 0 ; ixy < nxy ; ++ixy) imgInf->imgDiler[ixy + offset] = imgInf->img[ixy];
		/* Calculate dilation/erosion */
		imgInfDilation(&(imgInf->imgDiler[offset]), gptContext->dilerParmL[dilerN], gptContext);
		/* Canonical Definition */
		defcan(&(imgInf->imgDiler[offset]), &(imgInf->canDiler[offset]), gptContext->nxy);
		/* Make an edge direction image */
		roberts8(&(imgInf->imgDiler[offset]), &(imgInf->angDiler[offset]), gptContext->nx, gptContext->ny);
		/* Make template for distance of the same direction (if nndegdType >= 2) */
		mkNnDistTbl(&(imgInf->angDiler[offset]), &(imgInf->nndegdTblDiler[(dilerN - 1) * gptContext->nxyDir]), gptContext->nndegdContext);
		/* Make weight data */
		if (gptContext->evalType == EVALNNDEGDWEIGHT || gptContext->evalType == EVALNNDEGDSQRWEIGHT)
			calWeight(&(imgInf->imgDiler[offset]), &(imgInf->weightDiler[offset]), gptContext->nx, gptContext->ny, gptContext->workImgD1, gptContext->workImgD2, gptContext->workImgD3, gptContext->workImgD4);
	}
}

/* Dilate image for */
void imgInfDilation(int *img, int dilerPara, GptContext *gptContext) {
	int dilerType, ixy;
	double *workD1 = gptContext->workImgD1;
	double *workD2 = gptContext->workImgD2;
	double coef = DILCOEF;

	if (dilerPara == 0) return;
	else if (dilerPara > 0) dilerType = 4; /* Dilation for black background and Erosion for white background */
	else dilerType = 5;
	for (ixy = 0 ; ixy < gptContext->nxy ; ++ixy) workD1[ixy] = img[ixy];
	dilationIter(dilerType, workD1, workD2, gptContext->nx, gptContext->ny, coef, abs(dilerPara));
	for (ixy = 0 ; ixy < gptContext->nxy ; ++ixy) img[ixy] = (int) workD2[ixy];
}

/* Dilate image for */
void imgInfSmooth(int *img, int dilerPara, GptContext *gptContext) {
	int dilerType, dilerType2, ixy;
	double *workD1 = gptContext->workImgD1;
	double *workD2 = gptContext->workImgD2;
	double coef = DILCOEF;

	if (dilerPara == 0) return;
	else if (dilerPara > 0) {
		dilerType  = 4; /* Dilation for black background and Erosion for white background */
		dilerType2 = 5;
	}
	else {
		dilerType  = 5;
		dilerType2 = 4;
	}
	for (ixy = 0 ; ixy < gptContext->nxy ; ++ixy) workD1[ixy] = img[ixy];
	dilationIter(dilerType,  workD1, workD2, gptContext->nx, gptContext->ny, coef, abs(dilerPara));
	dilationIter(dilerType2, workD2, workD1, gptContext->nx, gptContext->ny, coef, abs(dilerPara));
	for (ixy = 0 ; ixy < gptContext->nxy ; ++ixy) img[ixy] = (int) workD1[ixy];
}

/* Dilate image for */
void imgInfSmoothD(int *img, double *imgSmooth, int dilerPara, double *workD1, int nx, int ny) {
	int dilerType, dilerType2, ixy, nxy = nx * ny;
	double coef = DILCOEF;

	if (dilerPara > 0) {
		dilerType  = 4; /* type = 4 : Dilation for black background and Erosion for white background */
		dilerType2 = 5;
	} else {
		dilerType  = 5;
		dilerType2 = 4;
	}
	for (ixy = 0 ; ixy < nxy ; ++ixy) imgSmooth[ixy] = img[ixy];

	if (dilerPara == 0) return;

	dilationIter(dilerType,  imgSmooth, workD1, nx, ny, coef, abs(dilerPara));
	dilationIter(dilerType2, workD1, imgSmooth, nx, ny, coef, abs(dilerPara));
}


/* Load image information for target image */
void loadImgInf(FILE *imgInfFp, ImgInf *imgInf, GptContext *gptContext) {
	size_t size;
	size = fread(imgInf->img,       sizeof(int),    gptContext->nxy,         imgInfFp);
	size = fread(imgInf->ang,       sizeof(int),    gptContext->nxy,         imgInfFp);
	size = fread(imgInf->can,       sizeof(double), gptContext->nxy,         imgInfFp);
	size = fread(imgInf->weight,    sizeof(double), gptContext->nxy,         imgInfFp);
#ifndef NOGPT
	size = fread(imgInf->gH0,       sizeof(double), gptContext->nxyDirNdGamma, imgInfFp);
	size = fread(imgInf->gH1x,      sizeof(double), gptContext->nxyDirNdGamma, imgInfFp);
	size = fread(imgInf->gH1y,      sizeof(double), gptContext->nxyDirNdGamma, imgInfFp);
#endif
	if (size == 0) {
		printf("Cannot read image information file at 0 in loadImgInf()\n");
		exit(0);
	}
	if (size == 0) {
		printf("Cannot read image information file at 1 in loadImgInf()\n");
		exit(0);
	}

	size = fread(imgInf->nndegdTbl, sizeof(double), gptContext->nxyDir,      imgInfFp);
	if (size == 0) {
		printf("Cannot read image information file at 2 in loadImgInf()\n");
		exit(0);
	}

	/* for dilation/Erosion data */
	if (gptContext->nDiler != 1) {
		size = fread(imgInf->imgDiler,       sizeof(int),    gptContext->nxy    * (gptContext->nDiler - 1), imgInfFp);
		size = fread(imgInf->angDiler,       sizeof(int),    gptContext->nxy    * (gptContext->nDiler - 1), imgInfFp);
		size = fread(imgInf->canDiler,       sizeof(double), gptContext->nxy    * (gptContext->nDiler - 1), imgInfFp);
		size = fread(imgInf->weightDiler,    sizeof(double), gptContext->nxy    * (gptContext->nDiler - 1), imgInfFp);
		size = fread(imgInf->nndegdTblDiler, sizeof(double), gptContext->nxyDir * (gptContext->nDiler - 1), imgInfFp);
	}
	if (size == 0) {
		printf("Cannot read image information file for diler in loadImgInf()\n");
		exit(0);
	}
}

/* Save image information for target image */
void saveImgInf(FILE *imgInfFp, ImgInf *imgInf, GptContext *gptContext) {
	fwrite(imgInf->img,       sizeof(int),    gptContext->nxy,         imgInfFp);
	fwrite(imgInf->ang,       sizeof(int),    gptContext->nxy,         imgInfFp);
	fwrite(imgInf->can,       sizeof(double), gptContext->nxy,         imgInfFp);
	fwrite(imgInf->weight,    sizeof(double), gptContext->nxy,         imgInfFp);
#ifndef NOGPT
	fwrite(imgInf->gH0,       sizeof(double), gptContext->nxyDirNdGamma, imgInfFp);
	fwrite(imgInf->gH1x,      sizeof(double), gptContext->nxyDirNdGamma, imgInfFp);
	fwrite(imgInf->gH1y,      sizeof(double), gptContext->nxyDirNdGamma, imgInfFp);
#endif
	fwrite(imgInf->nndegdTbl, sizeof(double), gptContext->nxyDir,      imgInfFp);

	/* for dilation/Erosion data */
	if (gptContext->nDiler != 1) {
		fwrite(imgInf->imgDiler,       sizeof(int),    gptContext->nxy    * (gptContext->nDiler - 1), imgInfFp);
		fwrite(imgInf->angDiler,       sizeof(int),    gptContext->nxy    * (gptContext->nDiler - 1), imgInfFp);
		fwrite(imgInf->canDiler,       sizeof(double), gptContext->nxy    * (gptContext->nDiler - 1), imgInfFp);
		fwrite(imgInf->weightDiler,    sizeof(double), gptContext->nxy    * (gptContext->nDiler - 1), imgInfFp);
		fwrite(imgInf->nndegdTblDiler, sizeof(double), gptContext->nxyDir * (gptContext->nDiler - 1), imgInfFp);
	}
}

/* Quick sort for integer array */
void quicksort_int(int *array, int *label, int left, int right) {
	int i, last;

	if (left >= right) return;
	swap_int(array, left, (left+right)/2);
	swap_int(label, left, (left+right)/2);
	last = left;
	for (i = left + 1; i <= right; ++i) {
		if (array[i] > array[left]) {
			++last;
			swap_int(array, last, i);
			swap_int(label, last, i);
		}
	}
	swap_int(array, left, last);
	swap_int(label, left, last);
	quicksort_int(array, label, left, last - 1);
	quicksort_int(array, label, last + 1, right);
}

/* Definite canonicalization */
void defcan(int *inImg, double *imgCan, int nxy) {
	int pos;
	double mean, var, varL1, ratio;

	calMeanVar(inImg, &mean, &var, &varL1, nxy);

	if (var == 0.0) var = 1.0;
	ratio = 1.0 / sqrt(var);
	for (pos = 0 ; pos < nxy ; ++pos) {
		imgCan[pos] = ratio * ((double) inImg[pos] - mean);
#ifdef NODEFCAN
		imgCan[pos] = inImg[pos];
#endif
#ifdef DEBUG
		printf("pos = (%d, %d) can = %f \n", pos % ROW, pos / ROW, imgCan[pos]);
#endif
	}
}

void calMeanVar(int *inImg, double *mean, double *var, double *varL1, int nxy) {
	int pos;
	double tmp;

	*mean = *var = *varL1 = 0.0;
	for (pos = 0 ; pos < nxy ; ++pos) *mean += (double) inImg[pos];
	*mean /= nxy;
	for (pos = 0 ; pos < nxy ; ++pos) {
		tmp  = (double) inImg[pos] - *mean;
		*var   += tmp * tmp;
		*varL1 += fabs(tmp);
	}
	return;
}

/* Extraction of gradient information by Roberts operator */
/* with 8-directional codes and strength */
void roberts8(int *img, int *imgAng, int nx, int ny) {
	static double ang5 =  7.0 / 8.0 * PI;
	static double ang4 =  5.0 / 8.0 * PI;
	static double ang3 =  3.0 / 8.0 * PI;
	static double ang2 =  1.0 / 8.0 * PI;
	static double ang1 = -1.0 / 8.0 * PI;
	static double ang0 = -3.0 / 8.0 * PI;
	static double ang7 = -5.0 / 8.0 * PI;
	static double ang6 = -7.0 / 8.0 * PI;
	int ix, iy, pos;  /* Pos variable */
	double delta_RD, delta_LD;
	double angle;

	for (pos = 0 ; pos < nx * ny ; ++pos) imgAng[pos] = NoDIRECTION;
	pos = 0;
	for (iy = 0 ; iy < ny - 1 ; ++iy) {
		for (ix = 0 ; ix < nx - 1 ; ++ix) {
			delta_RD = img[pos + 1] - img[pos + nx];
			delta_LD = img[pos] - img[pos + nx + 1];
			//if (delta_RD * delta_RD + delta_LD * delta_LD > 1.0) printf("direct strength = %f \n", sqrt(delta_RD * delta_RD + delta_LD * delta_LD));
			if (delta_RD * delta_RD + delta_LD * delta_LD > DIRECTEPS * DIRECTEPS) {
				if (fabs(delta_RD) < EPS) {
					if (delta_LD > 0) imgAng[pos] = 3;
					if (delta_LD < 0) imgAng[pos] = 7;
				} else {
					angle = atan2(delta_LD, delta_RD);
					if (angle > ang5) imgAng[pos] = 5;
					else if (angle > ang4) imgAng[pos] = 4;
					else if (angle > ang3) imgAng[pos] = 3;
					else if (angle > ang2) imgAng[pos] = 2;
					else if (angle > ang1) imgAng[pos] = 1;
					else if (angle > ang0) imgAng[pos] = 0;
					else if (angle > ang7) imgAng[pos] = 7;
					else if (angle > ang6) imgAng[pos] = 6;
					else imgAng[pos] = 5;
				}
			}
			++pos;
		}
		++pos; /* x < nx - 1 のため */
	}
}

/* Bluer by using extend image size (nx + 2) x (ny + 2) */
void bluer33Iter(double *extInImg, double *extWorkImg, int *pos33, double *w33, int nIter, int nx, int ny) {
	int loop, pos;

	for (loop = 0 ; loop < nIter ; loop += 2) {
		bluer33(extInImg, extWorkImg, pos33, w33, nx, ny);
		if (loop == nIter - 1) {
			for (pos = 0 ; pos < (nx + 2) * (ny + 2) ; ++pos) extInImg[pos] = extWorkImg[pos];
			break;
		}
		bluer33(extWorkImg, extInImg, pos33, w33, nx, ny);
	}
}

/* Bluer by using extend image size (nx + 2) x (ny + 2) */
void bluer33(double *extInImg, double *extOutImg, int *pos33, double *w33, int nx, int ny) {
	int ix, iy, ip, pos;
	double sum;
	for (iy = 0 ; iy < ny ; ++iy) {
		pos = (iy + 1) * (nx + 2) + 1;
		for (ix = 0 ; ix < nx ; ++ix) {
			sum = 0.0;
			for (ip = 0 ; ip < 9 ; ++ip) {
				sum += w33[ip] * extInImg[pos + pos33[ip]];
			}
			extOutImg[pos++] = sum;
		}
	}
}

/* Swap integers in int array */
void swap_int(int *array, int i, int j) {
	int temp;
	temp = array[i];
	array[i] = array[j];
	array[j] = temp;
}

/*
 * type = 0 : Gaussian
 *        1 : grad^2
 *        2 : dilation by max
 *        3 : erosion  by min
 *        4 : dilation by convergence (for ackground = 0)
 *        5 : erosion  by divergence
 *  Input image array is broken
 */
void dilationIter(int type, double *inImgD, double *outImgD, int nx, int ny, double coef, int nIter) {
	int ix, iy, jx, jy, ixy, iter;
	int nxy = nx * ny;
	double tmp, prop, w;
	// printf("type = %d nx = %d ny = %d coef = %f nIter %d \n", type, nx, ny, coef, nIter);
	if (nIter == 0) for (ixy = 0 ; ixy < nx * ny ; ++ixy) {
		outImgD[ixy] = inImgD[ixy];
		return;
	}
	if (type == 0) {
		/* Gaussian */
		for (iy = 0 ; iy < ny ; ++iy) {
			for (ix = 0 ; ix < nx ; ++ix) {
				prop = 0.0;
				tmp  = 0.0;
				for (jy = 0 ; jy < ny ; ++jy) {
					for (jx = 0 ; jx < nx ; ++jx) {
						prop += w = exp(-(coef * ((ix - jx) * (ix - jx) + (iy - jy) * (iy - jy))));
						tmp  += w * inImgD[jx + nx * jy];
					}
				}
				outImgD[ix + nx * iy] = tmp / prop;
				outImgD[ix + nx * iy] = tmp;
			}
		}
	} else if (type >= 2 && type <= 5){
		/* Iternation */
		for (iter = 0 ; iter < nIter ; iter += 2) {
			dilation(type, inImgD, outImgD,nx, ny, coef);
			if (iter == nIter - 1) break;
			dilation(type, outImgD, inImgD, nx, ny, coef);
		}
		if (nIter %2 == 0) 	for (ixy = 0 ; ixy < nxy ; ++ixy) outImgD[ixy] = inImgD[ixy];
	}
}

void dilation(int type, double *inImgD, double *outImgD, int nx, int ny, double coef) {
	int ix, iy, jx, jy, sx, sy, ex, ey;
	double  gradX, gradY, tmp;
	// printf("type = %d nx = %d ny = %d coef = %f nIter %d \n", type, nx, ny, coef, nIter);

	/* Dilation */

	for (iy = 0 ; iy < ny ; ++iy) {
		for (ix = 0 ; ix < nx ; ++ix) {
			sx = (ix - 1 >= 0)? (ix - 1) : 0; ex = (ix + 1 < nx)? (ix + 1) : (nx - 1);
			sy = (iy - 1 >= 0)? (iy - 1) : 0; ey = (iy + 1 < ny)? (iy + 1) : (ny - 1);
			tmp = inImgD[ix + nx * iy];
			switch (type) {
			case 1:
				gradX = (inImgD[ex + nx * iy] - inImgD[sx + nx * iy]) / (ex - sx);
				gradY = (inImgD[ix + nx * sy] - inImgD[ix + nx * ey]) / (ey - sy);
				tmp += coef * sqrt(gradX * gradX + gradY * gradY);
				break;
			case 2:
				/* Calculate maximum */
				for (jy = sy ; jy <= ey ; ++jy) {
					for (jx = sx ; jx <= ex ; ++jx) {
						tmp = (tmp >= inImgD[jx + nx * jy]) ? tmp : inImgD[jx + nx * jy];
					}
				}
				break;
			case 3:
				/* Calculate minimum */
				for (jy = sy ; jy <= ey ; ++jy) {
					for (jx = sx ; jx <= ex ; ++jx) {
						tmp = (tmp <= inImgD[jx + nx * jy]) ? tmp : inImgD[jx + nx * jy];
					}
				}
				break;
			case 4:
				gradX = hingeFunc(inImgD[sx + nx * iy] - inImgD[ix + nx * iy])
				      + hingeFunc(inImgD[ex + nx * iy] - inImgD[ix + nx * iy]);
				gradY = hingeFunc(inImgD[ix + nx * sy] - inImgD[ix + nx * iy])
                              		+ hingeFunc(inImgD[ix + nx * ey] - inImgD[ix + nx * iy]);
				//tmp += coef * (gradX + gradY);
				//tmp += coef * (gradX * gradX + gradY * gradY);
				//tmp += coef * sqrt(gradX * gradX + gradY * gradY);
				tmp += coef * (sqrt(gradX * gradX + gradY * gradY) + gradX + gradY);
				break;
			case 5:
				gradX = hingeFunc(inImgD[ix + nx * iy] - inImgD[sx + nx * iy])
	       			+ hingeFunc(inImgD[ix + nx * iy] - inImgD[ex + nx * iy]);
				gradY = hingeFunc(inImgD[ix + nx * iy] - inImgD[ix + nx * sy])
                              		+ hingeFunc(inImgD[ix + nx * iy] - inImgD[ix + nx * ey]);
				tmp -= coef * (sqrt(gradX * gradX + gradY * gradY) + gradX + gradY);
				break;
			}
			outImgD[ix + nx * iy] = tmp;
		}
	}
}

/* hingeFunc */
double hingeFunc(double in) {
	return (in >= 0) ? in : 0;
}

/* Copy an int image */
void imgCopy(int *img1, int *img2, int nxy) {
	int pos;
	for (pos = 0 ; pos < nxy ; ++pos) {
		img2[pos] = img1[pos];
	}
}

/* Copy a double image */
void canCopy(double *can1, double *can2, int nxy) {
	int pos;
	for (pos = 0 ; pos < nxy ; ++pos) can2[pos] = can1[pos];
}

/* Extend an int image */
void imgExtend(int *img1, int *img2, int nx, int ny, int nxExt, int nyExt) {
	int pos, ix, iy;
	for (pos = 0 ; pos < (nx + 2 * nxExt) * (ny + 2 * nyExt) ; ++pos) img2[pos] = 0;
	for (iy = 0 ; iy < ny ; ++iy) {
		pos = (iy + nyExt) * (nx + 2 * nxExt);
		for (ix = 0 ; ix < nx ; ++ix) {
			img2[pos] = img1[ix + nx * iy];
		}
	}
}

/* For debug */
/* Show a pair of intergers in an int array */
void imgPr(int *img1, int *img2, int nx, int ny) {
	int x, y;
	for (y = 0 ; y < ny ; ++y) {
		for (x  = 0 ; x < nx ; ++x) {
			printf("(%d,%d) img1 = %d, img2 = %d \n", x, y, img1[y * nx + x], img2[y * nx + x]);
		}
	}
}

/* Show GPT parameter */
void gptPr(double *gpt, char *st) {
	int i;
	printf("%s \n", st);
	for(i = 0 ; i < 3 ; ++i) {
		printf("%10.6f  %10.6f  %10.6f\n", gpt[i], gpt[i + 3], gpt[i + 6]);
	}
}

#ifdef DEBUGDILATE
int main(int argc, char *argv[]) {
	int nx = COL;
	int ny = ROW;
	int rgCatN   = 2;
	int rgNinCat = 3;
	GptContext *gptContext;
	char rgFile[128];
	char outFile[128];
	int nDil = 10, ixy;

	int img[nx * ny];
	double work1[nx * ny];
	double work2[nx * ny];

	gptContext  = (GptContext *) malloc(sizeof(GptContext));
	gptContext->nx = nx;
	gptContext->ny = ny;
	gptContext->nxy = nx * ny;
	gptContext->workImgD1 = work1;
	gptContext->workImgD2 = work2;


	sprintf(rgFile, "%s/%s/ldg%d_%04d_gray.pgm", WORKDIR, PGMDIR, rgCatN, rgNinCat);
	load_image_file(rgFile, img, nx, ny, 0);

	for (ixy = 0 ; ixy < nx * ny ; ++ixy) img[ixy] = 0;
	img[nx / 2 + nx * (ny / 2)] = WHITE;

	imgInfDilation(img, nDil, gptContext);
	sprintf(outFile, "%s/outImg/psfAv_%d_%04d_%1d.pgm", WORKDIR, rgCatN, rgNinCat, nDil);
	save_image_file(outFile, img, nx, ny);
	return 0;
}
#endif
