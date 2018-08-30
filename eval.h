/*
 * eval.h
 *
 *  Created on: 2013/07/02
 *      Author: yamasita
 */

#ifndef EVAL_H_
#define EVAL_H_

#define NNDEGDORG       0
#define NNDEGDTBL       1
#define NNDEGDTMPLT     2
#define NNDEGDTMPLTHALF 3

/* Element of scan table */
typedef struct {
  int    dx;     /* Displacement of x       */
  int    dy;     /* Displacement of y       */
  int    dPos;   /* Displacement index      */
  int    l1Dist; /* l1 distance from center */
  double dist;   /* l2 distance from center */
} ScanNode;

/* Dnn calculation context */
typedef struct {
  /* Fixed parameters */
  int         nx;            /* # of horizontal pixels */
  int         ny;            /* # of vertical pixels */
  int         nDir;          /* # of quantized directions  (= 8) */
  /* Values fixed by the above parameters */
  int         nst;           /* # of scan table */
  int         dDiagSq;       /* Squared length of diagonal of image */
  double      dDiag;         /* length of diagonal of image */
  /* Tables independent of content of image */
  int        *nndegdBound;   /* Distance to the boundary from each point in order to prune search */
  ScanNode   *scanNodes;     /* Scan table */
  /* Work region */
  /* Parameter for execution */
  int         nndegdType;     /* Type of calculation of the minimum distance of the same direction */
  /* Output*/
  double     nndegd;          /* The minimum distance of the same direction */
  /* Work Directory */
  char       *workDir;        /* Directory for work file */
} NndegdContext;

/* Prototype declaration */
double imgCor(double *img1, double *img2, int nxy);
void   calWeight(int *img, double *weight, int nx, int ny, double *workImgD1, double *workImgD2, double *workImgD3, double *workImgD4);
double diffImgEr(int *img1, int *img2, int nErosion, int nx, int ny, double *workImgD1, double *workImgD2, double *workImgD3, double *workImgD4, int *diffErImg);
void   mkScanTbl(NndegdContext *nndegdContext);
void   mkNnDistTbl(int *imgAng, double *nndegdTbl, NndegdContext *nndegdContext);
void   nndegdInit(int nx, int ny, NndegdContext *nndegdContext);
double nndegd(int type, int *g_ang1, int *g_ang2, NndegdContext *nndegdContext, double *nndegdTbl);
double nndegdEv(int type, int evType, int *ang1, int *ang2, double *weight1, double *weight2, NndegdContext *nndegdContext, double *nndegdTbl);
double wdchEv(int *ang1, int *ang2, int nx, int ny);
#endif /* EVAL_H_ */
