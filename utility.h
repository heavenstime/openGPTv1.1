/*
 * utility.h
 *
 *  Created on: 2012/11/02
 *      Author: yamasita
 */

/* Prototype declaration of functions */
void load_image_file(char *filename, int *img, int nx, int ny, int nExt);
void save_image_file(char *filename, int *img, int nx, int ny);

/* Load and save image information file */
void allcImgInf(ImgInf *imgInf, GptContext *gptContext);
void loadMkImgInf(char *fileName, ImgInf *imgInf, GptContext *gptContext);
void loadImgInf(FILE *imgInfFp, ImgInf *imgInf, GptContext *context);
void saveImgInf(FILE *imgInfFp, ImgInf *imgInf, GptContext *context);

/* Image dilation/erorsion */
void imgInfDilation(int *img, int dilerPara, GptContext *gptContext);
void imgInfSmooth(int *img, int dilerPara, GptContext *gptContext);
void imgInfSmoothD(int *img, double *imgSmooth, int dilerPara, double *workD1, int nx, int ny);

/* Image processing */
void roberts8(int *img, int *imgAng, int nx, int ny);
void defcan(int *inImg, double *imgCan, int nxy);
void calMeanVar(int *inImg, double *mean, double *var, double *varL1, int nxy);

/* Quick sort of integers */
void quicksort_int(int *array, int *label, int left, int right);

/* Bluer by 3 x 3 mask */
void bluer33Iter(double *extInImg, double *extWorkImg, int *pos33, double *w33, int nIter, int nx, int ny);

/* Dilation and erosion of image */
void dilationIter(int type, double *inImgD, double *outImgD, int nx, int ny, double coef, int nIter);

/* Copy image */
void imgCopy(int *img1, int *img2, int nxy);
void canCopy(double *can1, double *can2, int nxy);

/* Display transform coefficients  */
void gptPr(double *gpt, char *st);
/* Display values of pixels of image */
void imgPr(int *img1, int *img2, int nx, int ny);

