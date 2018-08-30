# OpenGPTv1

## Abstract
 * Numerical pattern recogition by GPT related method with kNN. It achieved 99.7% in the experiment with MNIST


## Features 
 * GAT or GPT can be selected.
 * Target image can be dilated/eroded and re-matched to transformed image.  
     (Mathematical morphology is used for dilation and erosion)
 *  Parameters can be changed in parameter.h.

## Directories 
 * `./gptWork/mnistPgm` : MNIST images are stored in pgm format.
 * `./gptWork/v1.1`     : template (by command = 5) and matching results are stored. 
 * `./gptWork/outImg`   : transformed images are stored for the image output command 

## File names in `./gptWork/v1.1/`
* (`execName` is a name to distinguish parameter setting, etc. For exmaple, mnistNNDEGD, mnistDIRECTEPS20, and so on)
* `execName_result_n_m`   : matching results between not more than 1,000 patters in trainind data and all test patterns. 
  * n and m are integers (For example `mnistNNDEGD_result_2_4`.)
  * n indicates the category of training data
  * m indicates the batch 
* `execName_rg_imgInf_n`  : Templates of training data of category n. (For example, `mnistNNDEGD_rg_imgInf_2`.)
* `execName_ts_imgInf_n`  : Templetes of test data of category n. (For example,  `mnistNNDEGD_ts_imgInf_2`.)

## How to use 
* `./OpenGPTv1 comamnd execName`
  *  comamnd =  0 : Match using all data  (Use after command 5) 
  *  command =  5 : Make image information file for registered and test images 
  *  command = 10 : Recognition by k-nearest neighbor and output error rate for k = 1, 2,...,10
  *  command = 11 : Show best 10 correlations of registered pattern for each test pattern 
  *  command = 12 : Show average of best 10 correlations of results 
  *  command = 14 : Make list of test pattern of which nearest pattern is in incorrect category 
* `./OpenGPTv1 comamnd  execName start_rg_digit end_rg_digit`
  *  comamnd =  1 : Match rg_digits from start_rg_digit to end_rg_digit to all ts_digits (Use after command 5) 
* `./OpenGPTv1 comamnd  execName start_rg_digit end_rg_digit rg_batch_start int rg_batch_end` 
  * comamnd =  1 : Match rg_digits from start_rg_digit to end_rg_digit to all ts_digits  (Use after command 5) 
* `./OpenGPTv1 comamnd execNme k`
  * command = 15 : Output ts image info for kNN error with top 10 rg image info (After command = 0)
* `./OpenGPTv1 comamnd` 
  * command = 20 : Match by a registration image and a test image (input images are specified by variable in program)
* `./OpenGPTv1 comamnd execName errFileName`
  * command = 21 : Match by a registration images and a test images described in list in file               
* `./OpenGPTv1 comamnd execName kNNErrorList` 
  * command = 22 : Output images of kNN error (Use with command 15)
* `./OpenGPTv1 comamnd execName`
  * command = 23 : Output dilated image (input images are specified by variable in program)

## Example 1 (Make tamplate, matching by a process, and show error rate of kNN)
 1. Parameter setting  
   edit `parameter.h`
 1. Bulid the program  
   `make`
 1. Make template  
   `./OpenGPTv1 5 mnistNNDEGD`
 1. Excute matching  
   `./OpenGPTv1 0 mnistNNDEGD`
 1. Excute kNN and output error rate for k = 1, 2,...,10  
   `./GPOpenGPTv1 10 mnistNNDEGD`

## Example 2 (Make tamplate by machine gpt00, matching by seprated 5 machines (gpt00, gpt01, gpt02, gpt03, gpt04), 
             and show error rate of kNN by machine gpt00. We assume all machines have the same gcc and liburaries)
1. Parameter setting at gpt00  
   edit `parameter.h`
1. Bulid the program at gpt00  
   `make`
1. Make template at gpt00  
   `./OpenGPTv1 5 mnistNNDEGD`
1. Distribute/copy files of gpt00 to gpt01, gpt02, gpt03, and gpt04 
1. Excute matching
   1. At gpt00  
     `./OpenGPTv1 1 mnistNNDEGD 0 1`
     1.  At gpt01  
     `./OpenGPTv1 1 mnistNNDEGD 2 3`
     1.  At gpt02  
     `./OpenGPTv1 1 mnistNNDEGD 4 5`
     1.  At gpt03  
     `./OpenGPTv1 1 mnistNNDEGD 6 7`
     1.  At gpt04  
     `./OpenGPTv1 1 mnistNNDEGD 8 9`
 1. Copy result files in ./gptWork/v1.1/ from gpt01, gpt02, gpt03, and gpt04 to gpt00
 1. Excute kNN and output error rate for k = 1, 2,...,10  
   `./OpenGPTv1 10 mnistNNDEGD`

## Example 3 (Make tamplate, matching by MPI with 5 cpu with 10 processes, and show error rate of kNN)
 *  This MPI program is only to submit 10 jobs for 10 category to 10 processes in order to reduce inputing command.
1. We prepare 6 comuters with NFS. 
   * One is a file server. Five are to calculate matching.
   * Home directory should be common.
   * OpenMPI has to be installed for each computer.
1. Describe host comfigulation file named hosts2-01-05.txt
   * The content of hosts2-01-05.txt  (2 processes runs for a computer)
---

 `gpt01.localdomain cpu=2`  
 `gpt02.localdomain cpu=2`  
 `gpt03.localdomain cpu=2`  
 `gpt04.localdomain cpu=2`  
 `gpt05.localdomain cpu=2`  

---

1. Parameter setting  
   `edit parameter.h`
 1.Bulid the program
   `make`  
   `make -f MakeMPI`
 1. Make template  
   `./OpenGPTv1 5 mnistNNDEGD`
 1. Excute matching by 6 computers (The following command is done at machine0)  
   `nohup mpirun -hostfile hosts2-01-05.txt -np 10 OpenGPTv1MPI mnistNNDEGD &`
 1. Excute kNN and output error rate for k = 1, 2,...,10  
   `./OpenGPTv1 10 mnistNNDEGD`

## Example 4 (After kNN, output every missclassified test image with its ten best matched training images)
 1. Make error list of kNN recognition (for example, let k = 8).  
   `./OpenGPTv1 15 mnistNNDEGD 8 > errorList.txt`
 1.  Output images to ./gptWork/outImg/  
   `./OpenGPTv1 22 mnistNNDEGD errorList.txt`

***********************************************************************************************
## Selections in parameter.h
* DIRMATCH    GPT/PPT integral with/without delta_{> 0}(q(nabla(f), q(nabla(g)))
  * If defined    :'with'  (This is standard)
  * If not define : without'

* GPTTYPE     Type of transformation
  * ALLGAT       : only GAT
  * ATLGATPPT    : GAT and PPT alternatively 
  * GATPPTAFTGAT : Only GAT for the first several iterations, and GAT and PPT alternatively/

* GPTDIR      Images on which GAT/PPT is applied. 
  * GPTFORTEST    : Only test patterns are transformed.         
  * GPTFORTESTREG : Test and training patters are transformed.

* GPTMETHOD   Method to calculate GPT/PPT
  * GPTORG        : Original method. (Computational complexity is the squre of number of pixeds)
  * GPTTMPLT      : Template method. (Computational complexity is the number of pixeds)

* GAUSSMETHOD Calculation of Gaussian function
  * GAUSSTRUNC   : Truncate calculation of Gaussian function
  * NOGUASSTRUNC : Do not truncate the  calculation of Gaussian function

* NNDEGDTYPE   Calculation type of NNDEGD to fix the parameter of Gaussian fuction.
  * NNDEGDORG       : Original calculation method (Computational complexity is the squre of number of pixeds)
  * NNDEGDTBL       : Seach from nearer points
  * NNDEGDTMPLT     : Template method (Computational complexity is the number of pixeds)
  * NNDEGDTMPLTHALF : Template method, only distance from test image to training image is caluclated 

* EVALTYPE     Matching measure for recognition
  * EVALCOR          : Canonical correlation
  * EVALNNDEGD       : NNDEGD
  * EVALNNDEGDVAR    : NNDEGD variance 
  * EVALNNDEGDWEIGHT : NNDEGD with curvature weight
  * EVALNNDEGDSQR    : Squared NNDEGD
  * EVALDIFFER       : Difference
  * EVALWDCH         : Histogram of edge direction with Gassian weight

* EVALNNDEGDTYPE    Calculation type of NNDEGD for matching measure
  * NNDEGDORG       : Original calculation method (
  * NNDEGDTBL       : Seach from nearer points
  * NNDEGDTMPLT     : Template method 
  * NNDEGDTMPLTHALF : Template method, only distance from test image to training image is caluclated 

## Selections in parameter.h

* TIMESFSTGAT   : number of first GAT calculation whenr GPTTYPE is set to GATPPTAFTGAT
* MAX_ITER      : number of iterations for GPT/PPT calc.
* STOP_ITER     : If corNew < STOP_ITER * col, stop loop

* GAMMALIST     : List of values for Gaussian fuction parameter when template method of GAT/PPT is chosen. 
* NGAMMA        : number of parameters in the list

* NDIFFER       : number of iteration of erosion for eroded difference evaluation 

* NDILER        : The number of types of dilation/erosion for training patters. (= 1 : Dilation/erosion is not applied)
* DILER         : List of numbers of iterations for dilation/erosion process 
* DILCOEF       : Step size for differential equation of dilation/erosion 
* DILCUTTHRESH  : Threshold to apply matching with dilated/eroded patterns. They are used when eval > DILCUTTHRESH 

* CURVSMOOTH    : number of iterations to smooth the image to calcuate curvature for weighed NNDEGD
* CURVITER      : number of iterations to smooth the covariance matrix component
* CURVCOEFF     : Weight for curvature in curvature weighed NNDEGD


* INVERSECHK    : Inverse check to exculde the case when an image are truncated at boundary but correlation is high
  * 0 : Skip inverse check
  * 1 : Do inverse check.
* SELFCORLIMIT  : Threshold of correlation of original and inverse transformed after transformed image 

* TRANSCHK      : Chech transformation to excude the case the transformation is too large.
  * 0 : Skip transform check
  * 1 : Do transform check */
* TRANSLIMIT    : Lower bound of det(A A^T) / tr(A A^T)/

* NEXT          : number of pixels for extension of peripheral region input image extension
* INPUTSMOOTH   : Smoothing parameter for all input data (number of iterations)
* WGT           : Coefficient for Gauss integral parameter multiplied to NNDEGD
* EPS           : Value assumed to be zero
* DIRECTEPS     : Value assumed to be zero for robert8 edge detector
* BLUERCOEFLIST : Coefficient list for 3x3 bluer coefficients to calculate curvature
* WDCHGAUSS     : Gaussian fuction of 5 x 5 region for evaluation with histogram of edge direction with Gassian weight
* WDCHDIRSMOOTH : Parameters for smoothing direction for evaluation with histogram of edge direction with Gassian weight

## Matlab/Octave program to convert original MNIST data to pgm files.
* `convData.m`    : Program to convert original MNIST data to pgm files. 
  * Set variables to specify directories for orinal data and output data.
* `writePgm.m`    : Used in convData.m to output a pgm file