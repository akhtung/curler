//The Expectation Maximum Algorithm
//the plot is right 
//add the weight to the optics part
//the microcluster is a core microcluster when the sum of its neighbors' weight exceed the minimum density
//add in competing learning
//the graph seems to be right, but not good enough
//try to reduce the usage of memory
//different out file

#include <fstream.h>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>

#include "mytime.h"

#define FALSE 0
#define TRUE 1
#define eigeneps 0.000001
#define pagesize 50
#define minprob 0.001
#define compressratio 1
#define compressratio2 1
#define iterationtime 5
#define BIGENOUGH_FLOAT	10000

//iterationtimeSetting: quadratic(3),cubic(3),ovalhelixes(10),iris(8),image(8),aml(3),Iyer(20)

struct S_BigEnoughItem
{
	int	ItemNO;
	short *ItemList;
};

class EM{
public:
	EM(char *ufile, int pnumcluster, int pdimension, int ppagecap);
	~EM();
	int pagecap;	
	void EMCluster(char *pinfile, char *poutfile, char *porderfile, float pepsilon, float pdistepsilon, float pmindensity);
private:
	int processedcluster;
	float * eigmem;
	float *PRlx;
	float *uXXT; 
	float *uX;
	float * tempmem;
	float *resultmem;
	float *N;
	float Det(float *matrix, int pdim);
	float setcoredistance(float pN, float pmindensity, int pneighcount, float *pneighdif, int *pneighindex);
	void popseedlist(int *pseedlist, float *pseeddif, int *pcluster2seed, int seedcount);
	void expandclusterorder(int pclusternum, float pdisteps, float pmindenstiy, int *pprocessed, FILE *pofp);
	void orderseedlist(int pseedcount, int *pseedlist, float *pseeddif, int *pcluster2seed);
	int updatedif(int *pseedflag, int *pcluster2seed, int *pseedlist, float *pseeddif, int *pneighindex, float *pneighdif, int pseedcount, int pneighcount);
	int neighbour(int clusternum, float pdistepsi, int *pneighindex, float *pneighdif, int *pprocessed);
	void optics(char *poutfile, float pdistepsi, float pmindensity);
	float probdifference(int numcluster1, int numcluster2);
	void filewrite(char *poutfile, float *pu, float *psigma, int pdim, int pnumcluster);
    int dim;
	int numcluster;
	float *RS; 
	char *infile, *outfile;
	float weight;
	float *u;
	float *sigma;
	float *PRxl;
	float *determinant;
	char * clusterfile;
	int sigmasize;
	float *dist;
	int samplesize;
	S_BigEnoughItem *tmp_ItemList;

	void eigenvalue(FILE *pofp, float *matrix, int pdim);
    void InverseMatrix(float *pmatrix,int pmatrixdim, int pnumcluster);
};

void main(int argc, char **argv) {
	char *minfile, *moutfile, *meanfile, *orderfile;
	int mdimension;
	int mnumcluster;
	int mpagecap;
	float mepsilon, mdistepsilon, mmindensity;
	Timer excutetime;
	
	printf("\nClustering with EM algorithm\n");
	/*
	if(argc != 11) { 
		printf("  ** Usage: %s <in-file> <num of cluster> <dimension> <page capacity> <epsilon> <minimum density> <distepsilon> <meanfile> <out-file> <ordered file>\n\n",argv[0]); exit(-2); 
	}
	minfile = argv[1];
//printf("the input file: %s \n", minfile);
	mnumcluster = atoi(argv[2]);
//	printf("the number of clusters: %d \n",mnumcluster);
	mdimension = atoi(argv[3]);
//	printf("the number of dimesion: %d \n",mdimension);
	mpagecap = atoi(argv[4]);
//	printf("the pagecap of the file: %d \n",mpagecap);
	mepsilon = atof(argv[5]);
//	printf("the epsilon: %f \n",mepsilon);
	mmindensity = atof(argv[6]);
//	printf("the mininum density: %f \n", mmindensity);
	mdistepsilon = atof(argv[7]);
//	printf("the distepsilon: %f \n", mdistepsilon);
	meanfile = argv[8];
//	printf("the mean file: %s\n",meanfile);
	moutfile = argv[9];
//	printf("the output file: %s\n",moutfile);
	orderfile = argv[10];
//	printf("the ordered file: %s\n",orderfile);
*/

/*
    minfile="Release\\image\\image.ascii";
	mnumcluster=500;
	mdimension=19;
	mpagecap=2310;
	//mp=0.5;
	mmindensity = 0;
	mepsilon=10;
	mdistepsilon = 0;
	moutfile="Release\\image\\image.out";
	meanfile="Release\\image\\image_500.ascii";
	orderfile="Release\\image\\image_500.clu"; */

    /* quadratic sample 
	minfile="sample.ascii";
	mnumcluster=186;
	mdimension=2;
	mpagecap=3720;
	mepsilon=0;
	mmindensity = 20;
	mdistepsilon = 0;
	meanfile="sample_seed.ascii";
	moutfile="sample_out.ascii";
	orderfile="sample_clu.ascii";*/

	/* quadratic 
	minfile="quadratic.txt";
	mnumcluster=150;
	mdimension=2;
	mpagecap=2354;
	mepsilon=0;
	mmindensity = 0;
	mdistepsilon = 0;
	meanfile="quadratic_seed.txt";
	moutfile="quadratic_out.ascii";
	orderfile="quadratic_clu.ascii";*/

    /* cubic 
	minfile="cubic589.ascii";
	mnumcluster=150;
	mdimension=2;
	mpagecap=589;
	mepsilon=0;
	mmindensity = 20;
	mdistepsilon = 0;
	meanfile="cubic_seed.txt";
	moutfile="cubic_out.ascii";
	orderfile="cubic_clu.ascii";*/

	/* sin 
	minfile="sin.txt";
	mnumcluster=250;
	mdimension=2;
	mpagecap=2500;
	mepsilon=0;
	mmindensity = 0;
	mdistepsilon = 0;
	meanfile="sin_seed.txt";
	moutfile="sin_out.ascii";
	orderfile="sin_clu.ascii";*/

    /* helix  
	minfile="ovalhelixes.txt";
	mnumcluster=450;
	mdimension=9;
	mpagecap=3000;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="ovalhelixes_seed450_ordered.txt";
	moutfile="ovalhelixes_out.ascii";
	orderfile="ovalhelixes_clu.ascii"; */

	/* iris 
	minfile="iris.ascii";
	mnumcluster=150;
	mdimension=4;
	mpagecap=150;
	mepsilon=0;
	mmindensity = 50;
	mdistepsilon = 0;
	meanfile="iris_seed150.ascii";
	moutfile="iris_out.ascii";
	orderfile="iris_clu.ascii";*/

	/* image 
	minfile="imagenew.ascii";
	mnumcluster=500;
	mdimension=16;
	mpagecap=2310;
	mepsilon=0;
	mmindensity = 30;
	mdistepsilon = 0;
	meanfile="imagenew_seed500.ascii";
	moutfile="imagenew_out.ascii";
	orderfile="imagenew_clu.ascii";*/

    /* combined
    minfile="combined4628.ascii";
	mnumcluster=594;
	mdimension=3;
	mpagecap=4628;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="combined594.ascii";
	moutfile="combined_out.ascii";
	orderfile="combined_clu.ascii";*/

    /* compresscombined3d 
    minfile="compresscombined3d.txt";
	mnumcluster=500;
	mdimension=3;
	mpagecap=3628;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="compresscombined3d_seed500.txt";
	moutfile="compresscombined3d_out.ascii";
	orderfile="compresscombined3d_clu.ascii";*/

    /* compresscombined4d (seed800 and seed200) 
    minfile="compresscombined4d.txt";
	mnumcluster=200;
	mdimension=4;
	mpagecap=3628;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="compresscombined4d_seed200.txt";
	moutfile="compresscombined4d_out.ascii";
	orderfile="compresscombined4d_clu.ascii";*/

    /* ovalcombined3d seed800 seed200 seed500 
    minfile="ovalcombined3d.txt";
	mnumcluster=800;
	mdimension=3;
	mpagecap=3000;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="ovalcombined3d_seed800.txt";
	moutfile="ovalcombined3d_out.ascii";
	orderfile="ovalcombined3d_clu.ascii";*/

    /* ovalcombined4d seed800 seed200 seed500 
    minfile="ovalcombined4d.txt";
	mnumcluster=500;
	mdimension=4;
	mpagecap=3000;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="ovalcombined4d_seed500.txt";
	moutfile="ovalcombined4d_out.ascii";
	orderfile="ovalcombined4d_clu.ascii"; */

    /* noisecombined3d seed500 
    minfile="noisecombined3d.txt";
	mnumcluster=500;
	mdimension=3;
	mpagecap=3000;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="noisecombined3d_seed500.txt";
	moutfile="noisecombined3d_out.ascii";
	orderfile="noisecombined3d_clu.ascii"; */

    /* combined7dnoise05 seed500 seed400 seed300 
    minfile="combined7dnoise05.txt";
	mnumcluster=300;
	mdimension=7;
	mpagecap=3000;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="combined7dnoise05_seed300.txt";
	moutfile="combined7dnoise05_out.ascii";
	orderfile="combined7dnoise05_clu.ascii"; */

    /* combined7dnoise30 seed500 seed400 seed300 
    minfile="combined7dnoise30.txt";
	mnumcluster=500;
	mdimension=7;
	mpagecap=3000;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="combined7dnoise30_seed500.txt";
	moutfile="combined7dnoise30_out.ascii";
	orderfile="combined7dnoise30_clu.ascii"; */

    /* testcombined6d seed500 seed400 seed300 
    minfile="testcombined6d.txt";
	mnumcluster=500;
	mdimension=6;
	mpagecap=3000;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="testcombined6d_seed500.txt";
	moutfile="testcombined6d_out.ascii";
	orderfile="testcombined6d_clu.ascii"; */

    /* testcombined5d seed500 seed400 seed300 
    minfile="testcombined5d.txt";
	mnumcluster=500;
	mdimension=5;
	mpagecap=3000;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="testcombined5d_seed500.txt";
	moutfile="testcombined5d_out.ascii";
	orderfile="testcombined5d_clu.ascii"; */

    /*
    minfile="aml72grainratiogene10n5.ascii";
	mnumcluster=72;
	mdimension=10;
	mpagecap=72;
	mepsilon=0;
	mmindensity = 0;
	mdistepsilon = 0;
	meanfile="aml72grainratiogene10n5_seed72.ascii";
	moutfile="aml_out.ascii";
	orderfile="aml_clu.ascii";*/


    /* Iyer517V2_seed360 
    minfile="Iyer517V2.ascii";
	mnumcluster=360;
	mdimension=18;
	mpagecap=517;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="Iyer517V2_seed360.ascii";
	moutfile="Iyer_out.ascii";
	orderfile="Iyer_clu.ascii";*/
	
	/* wdbc */
    minfile="wdbcn5.txt";
	mnumcluster=200;
	mdimension=30;
	mpagecap=569;
	mepsilon=0;
	mmindensity = 10;
	mdistepsilon = 0;
	meanfile="wdbcn5_seed200.txt";
	moutfile="wdbc_out.txt";
	orderfile="wdbc_clu.txt";

	excutetime.start();
	EM em(meanfile, mnumcluster, mdimension, mpagecap);
	em.EMCluster(minfile, moutfile, orderfile, mepsilon, mdistepsilon, mmindensity);	
	excutetime.end();
}

EM::EM(char *ufile, int pnumcluster, int pdimension, int ppagecap){
	int i=0,j=0,k=0;
	FILE *ufp;
	float meanval;
	int index;
	
	/* set the members of EM */

	sigmasize = pdimension * (pdimension + 1) / 2;
	samplesize = pnumcluster * compressratio;//*pnumcluster;

	pagecap = ppagecap;
	dim = pdimension;
	numcluster = pnumcluster;
	clusterfile = "sample1_170.clu";
		
	RS = (float *) malloc (sizeof (float) * dim);
	u = (float * ) malloc ( sizeof (float) * pnumcluster * dim);
	sigma = (float * ) malloc (sizeof (float) * pnumcluster * sigmasize);
	determinant = (float * ) malloc (sizeof (float) * pnumcluster);
	N = (float *) malloc (sizeof(float) * pnumcluster);
	tempmem = (float *) malloc (sizeof(float) * dim * dim);
	resultmem = (float *) malloc (sizeof(float) * dim * dim);
	eigmem = (float *) malloc (sizeof(float) * dim * (dim+4));

	//test the memory allocation
	if (RS == NULL){
		printf("Not enough memory for RS!!\n");
		exit(0);
	}
	if (u == NULL){
		printf("Not enough memory for u!!\n");
		exit(0);
	}
	if (sigma == NULL){
		printf("Not enough memory for sigma!!\n");
		exit(0);
	}
	if (determinant == NULL){
		printf("Not enough memory for determinant!!\n");
		exit(0);
	}
	if (N == NULL){
		printf("Not enough memory for N!!\n");
		exit(0);
	}
	if (tempmem == NULL){
		printf("Not enough memory for tempmem!!\n");
		exit(0);
	}
	if (eigmem == NULL){
		printf("Not enough memory for eigmem!!\n");
		exit(0);
	}

	/* open the mean value file */
	if ((ufp = fopen(ufile, "r"))==NULL)
	{
		printf("Mean File Not Opened!");
		exit(0);
	}

	/* read the initial mean value */
	for (i = 0; i < pnumcluster; i++)
	{
		for (j = 0; j < dim; j++){
			fscanf(ufp, "%f\t", &meanval);
			u[i * dim + j] = meanval;
		}
		fscanf(ufp, "\n");
	}
		
	/* initialize the sigma matrix */
	
	for (i = 0; i < pnumcluster; i++){	
		index = 0;
		for (j = 0; j < dim; j++){
			for (k = 0; k < j + 1; k++){				
				if (k == j) {
					sigma[i * sigmasize + index] = 1;				
					index = index + 1;
				}
				else {
					sigma[i * sigmasize + index] = 0;
					index = index + 1;
				}
			}
		}

	}
}

EM::~EM(){
    //free(RS);
	//free(u);
	//free(sigma);
	//free(determinant);
	free(N);
	free(tempmem);
	free(resultmem);
	free(eigmem);
}

void EM::InverseMatrix(float *pmatrix, int pmatrixdim, int pnumcluster){
	float *tempmatrix;
	int i,j,k,m;
	int count;
	float temp;
	int mat_index;
	int matrixsize;
	float *resultmatrix;
	short *IS, *JS;
	double largest;

	matrixsize = pmatrixdim * (pmatrixdim + 1) / 2;
	/*0                  0*1/2+0 */
	/*1, 2               1*2/2+0 or 1 */
	/*3, 4, 5            2*3/2+0 or 1 or 2 */
	/*6, 7, 8, 9         3*4/2+0 or 1 or 2 or 3 */
	count = pmatrixdim * pmatrixdim * pnumcluster;

	tempmatrix = tempmem;
	resultmatrix = resultmem;
//	tempmatrix = (float *) malloc (sizeof(float) * pmatrixdim * pmatrixdim);
//	resultmatrix = (float *) malloc (sizeof(float) * pmatrixdim * pmatrixdim);
	IS = (short *) malloc (sizeof(short) * pmatrixdim);
	JS = (short *) malloc (sizeof(short) * pmatrixdim);

	/* the number of clusters stored in pmatrix is pnumcluster */
	for (i = 0; i < pnumcluster; i++){
		for (j = 0; j < pmatrixdim; j++){
			IS[j] = j;
			JS[j] = j;
		}

		/* tempmatrix: read the data for the current cluster from pmatrix */
		for (j = 0; j < pmatrixdim; j++)
			for (k = 0; k < pmatrixdim; k++){
				if (j < k)
					mat_index = (k + 1) * k / 2 + j;
				else 
					mat_index = (j + 1) * j / 2 + k;
				tempmatrix[j * pmatrixdim + k] = pmatrix[i * matrixsize + mat_index];
			}

	    /* resultmatrix: set the element at diagonal as 1, else as 0 */
		for (j = 0; j < pmatrixdim; j++)
			for (k = 0; k < pmatrixdim; k++)
				if (j == k)
					resultmatrix[j * pmatrixdim + k] = 1;
				else
					resultmatrix[j * pmatrixdim + k] = 0;

		/* tempmatrix: move the largest element to the diagonal */
		for (j = 0; j < pmatrixdim; j++){

			/* find the largest element in tempmatrix */
			largest = fabs(tempmatrix[IS[j] * pmatrixdim + JS[j]]);
			for (k = 0; k < pmatrixdim; k++)
				for (m = 0; m < pmatrixdim; m++)
					if (fabs(tempmatrix[k * pmatrixdim + m]) > largest){
						largest = fabs(tempmatrix[k * pmatrixdim + m]);
						IS[j] = k;
						JS[j] = m;
					}

			/* exchange the rows */
			if (IS[j] != j){
				for (k = 0; k < pmatrixdim; k++){
					temp = tempmatrix[IS[j] * pmatrixdim + k];
					tempmatrix[IS[j] * pmatrixdim + k] = tempmatrix[j * pmatrixdim + k]; 
					tempmatrix[j * pmatrixdim + k] = temp;
				}

			}

			/* exchange the cols */
			if (JS[j] != j){
				for (k = 0; k < pmatrixdim; k++){
					temp = tempmatrix[JS[j] + k * pmatrixdim];
					tempmatrix[JS[j] + k * pmatrixdim] = tempmatrix[j + k * pmatrixdim]; 
					tempmatrix[j  + k * pmatrixdim] = temp;
				}
			}

            /* normalize row j */
			tempmatrix[j * pmatrixdim + j] = 1 / tempmatrix[j * pmatrixdim + j];
			for (k = 0; k < pmatrixdim; k++)
				if (k != j)
					tempmatrix[j * pmatrixdim + k] = tempmatrix[j * pmatrixdim + k] * tempmatrix[j * pmatrixdim + j];
			
			/* adjust elements not in row j and column j */
			for (k = 0; k < pmatrixdim; k++)
				for (m = 0; m < pmatrixdim; m++)
					if ((k != j) && (m != j))
						tempmatrix[k * pmatrixdim + m] = tempmatrix[k * pmatrixdim + m] -
							tempmatrix[k * pmatrixdim + j] * tempmatrix[j * pmatrixdim + m];

		    /* adjust elements in column j */
			for (k = 0; k < pmatrixdim; k++)
				if (k != j)
					tempmatrix[k * pmatrixdim + j] = - tempmatrix[k * pmatrixdim + j] * tempmatrix[j * pmatrixdim + j];		
		}//j<pmatrixdim

		for (j = pmatrixdim-1; j >= 0;j--){
			/* exchange the rows */
			if (JS[j] != j){
				for (k = 0; k < pmatrixdim; k++){
					temp = tempmatrix[JS[j] * pmatrixdim + k];
					tempmatrix[JS[j] * pmatrixdim + k] = tempmatrix[j * pmatrixdim + k]; 
					tempmatrix[j * pmatrixdim + k] = temp;
				}

			}

			/* exchange the cols */
			if (IS[j] != j){
				for (k = 0; k < pmatrixdim; k++){
					temp = tempmatrix[IS[j] + k * pmatrixdim];
					tempmatrix[IS[j] + k * pmatrixdim] = tempmatrix[j + k * pmatrixdim]; 
					tempmatrix[j  + k * pmatrixdim] = temp;
				}
			}
		}// j>=0

        /* resultmatrix=>pmatrix */
		for (j = 0; j < pmatrixdim; j++)
			for (k = 0; k <=j; k++)
				pmatrix[i * matrixsize + j * (j + 1) / 2 + k] = resultmatrix[j * pmatrixdim + k];

	}//for i<pnumcluster

//	free(tempmatrix);
//	free(resultmatrix);
	free(IS);
	free(JS);
}

void EM::eigenvalue(FILE *pofp, float *matrix, int pdim){
	int i, j, k;
	double maxelement;
	int p, q;
	float *S;
	double thita;
	int sign;
	float temp;
	int matrix_index, matrix_index1, matrix_index2;
		
	float m,n;
	float *temprowp, *temprowq, *tempcolp, *tempcolq;

	S = eigmem;
	temprowp = &eigmem[pdim * pdim];
	temprowq = &eigmem[pdim * pdim + pdim];
	tempcolp = &eigmem[pdim * pdim + pdim * 2];
	tempcolq = &eigmem[pdim * pdim + pdim * 3];

	/* step 1:set S to be the identity matrix; */
	for (i=0; i<pdim; i++)
		for (j=0; j<pdim; j++)
			if (i==j) S[i*pdim+j]=1;
			else S[i*pdim+j]=0;

	/* step2: get the maximun non-diagonal element apq of the matrix; */

	/* 00 */
    /* 10 11 */
	/* 20 21 22 */
	/* 30 31 32 33 */
	/* ... */
	maxelement=fabs(matrix[1]);
	p=0; /* column index */
	q=1; /* row index */

	for (i=2; i<pdim; i++)
		for (j=0; j<i; j++){
			matrix_index = i * (i + 1) / 2 + j;

			if (fabs(matrix[matrix_index])>maxelement){
				p=j;
				q=i;
				maxelement=fabs(matrix[matrix_index]);
			}
		}

	/* step3:if |apq|<eigeneps, terminates; */

	/* step4: calculate sin(2*thita), sin(thita), cos(thita); */
	
	while (fabs(maxelement)>eigeneps){
		if (p > q){//?, we know p<q since j<i
			matrix_index = p * (p + 1) / 2 + q;
		}
		else{
			matrix_index = q * (q + 1) / 2 + p;
		}
		m=-matrix[matrix_index];
		matrix_index1 = p * (p + 1) / 2 + p;/*diagonal element with row index p and col index p*/
		matrix_index2 = q * (q + 1) / 2 + q;/*diagonal element with row index q and col index q*/
		n=(matrix[matrix_index2]-matrix[matrix_index1])/2;
		if (n>=0) sign=1;
		else sign=-1;
		thita=asin(sign*m/sqrt(pow(m,2)+pow(n,2))) / 2;
			
		//step5: calculate the new element for the matrix;
		for (i=0; i<pdim; i++){
			//bpp,bpq,bqp
			if (i==p){
				matrix_index1 = p * (p + 1) / 2 + p;
				matrix_index2 = q * (q + 1) / 2 + q;
				if (p < q)
					matrix_index = (q + 1) * q / 2+ p;
				else 
					matrix_index = (p + 1) * p / 2+ q;
				temprowp[i]=float(matrix[matrix_index1]*pow(cos(thita),2)
					+matrix[matrix_index2]*pow(sin(thita),2)+matrix[matrix_index]*sin(2*thita));
				tempcolp[i]=temprowp[i];
				temprowq[i]=float((matrix[matrix_index2]-matrix[matrix_index1])*sin(2*thita)/2+matrix[matrix_index]*cos(2*thita));
				tempcolq[i]=temprowq[i];				
			}

			//bqq,bpq,bqp
			else if (i==q){
				matrix_index1 = p * (p + 1) / 2 + p;
				matrix_index2 = q * (q + 1) / 2 + q;
				if (p < q)
					matrix_index = (q + 1) * q / 2+ p;
				else 
					matrix_index = (p + 1) * p / 2+ q;
				temprowq[i]=float(matrix[matrix_index1]*pow(sin(thita),2)
					+matrix[matrix_index2]*pow(cos(thita),2)-matrix[matrix_index]*sin(2*thita));
				tempcolq[i]=temprowq[i];
				temprowp[i]=float((matrix[matrix_index2]-matrix[matrix_index1])*sin(2*thita)/2+matrix[matrix_index]*cos(2*thita));
				tempcolp[i]=temprowp[i];					
			}

			else{
				if (p > i)
					matrix_index1 = p * (p + 1) / 2 + i;
				else
					matrix_index1 = i * (i + 1) / 2 + p;
				if (q > i)
					matrix_index2 = q * (q + 1) / 2 + i;
				else
					matrix_index2 = i * (i + 1) / 2 + q;

				temprowp[i]=float(matrix[matrix_index1]*cos(thita)+matrix[matrix_index2]*sin(thita));
				temprowq[i]=float(-matrix[matrix_index1]*sin(thita)+matrix[matrix_index2]*cos(thita));

				if (p > i)
					matrix_index1 = p * (p + 1) / 2 + i;
				else
					matrix_index1 = i * (i + 1) / 2 + p;

				if (q > i)
					matrix_index2 = q * (q + 1) / 2 + i;
				else
					matrix_index2 = i * (i + 1) / 2 + q;

				tempcolp[i]=float(matrix[matrix_index1]*cos(thita)+matrix[matrix_index2]*sin(thita));
				tempcolq[i]=float(-matrix[matrix_index1]*sin(thita)+matrix[matrix_index2]*cos(thita));
			}
		}//for i<pdim

		//update matrix A;
		for (i=0; i<pdim; i++){
			if (p > i)
				matrix_index1 = p * (p + 1) / 2 + i;
			else
				matrix_index1 = i * (i + 1) / 2 + p;

			if (q > i)
				matrix_index2 = q * (q + 1) / 2 + i;
			else
				matrix_index2 = i * (i + 1) / 2 + q;

			matrix[matrix_index1]=temprowp[i];
			matrix[matrix_index2]=temprowq[i];

			if (p > i)
				matrix_index1 = p * (p + 1) / 2 + i;
			else
				matrix_index1 = i * (i + 1) / 2 + p;

			if (q > i)
				matrix_index2 = q * (q + 1) / 2 + i;
			else
				matrix_index2 = i * (i + 1) / 2 + q;

			matrix[matrix_index1]=tempcolp[i];
			matrix[matrix_index2]=tempcolq[i];
		}

		/* step6: update S=S*R(p,q,thita) */
		for (i=0; i<pdim; i++){
			tempcolp[i]=float(S[i*pdim+p]*cos(thita)+S[i*pdim+q]*sin(thita));
			tempcolq[i]=float(S[i*pdim+p]*(-sin(thita))+S[i*pdim+q]*cos(thita));
		}

		for (i=0; i<pdim; i++){
			S[i*pdim+p]=tempcolp[i];
			S[i*pdim+q]=tempcolq[i];
		}
		
		/* step2: get the maximun non-diagonal element apq of the matrix; */
		maxelement=fabs(matrix[1]);
		p=0;
		q=1;
	
		for (i=0; i<pdim; i++)
			for (j=0; j<i; j++){
				matrix_index = i * (i + 1) / 2 + j;
				if (fabs(matrix[matrix_index])>maxelement){
					p=j;
					q=i;
					maxelement=fabs(matrix[matrix_index]);
				}
			}
	}//while (fabs(maxelement)>eigeneps)		

	/* sort the eigenvalues */
	for(i = 0; i < pdim; i++)
		for (j = 0; j < pdim - i -1; j++)			
			if (fabs(matrix[j * (j + 1) / 2 + j]) > fabs(matrix[(j + 1) * (j + 2) / 2 + (j + 1)])){
				temp = matrix[j * (j + 1) / 2 + j];
				matrix[j * (j + 1) / 2 + j] = matrix[(j + 1) * (j + 2) / 2 + (j + 1)];
				matrix[(j + 1) * (j + 2) / 2 + (j + 1)] = temp;
				for (k = 0; k < pdim; k++){
					temp = S[k * pdim + j];
					S[k * pdim + j] = S[k * pdim + j + 1];
					S[k * pdim + j + 1] = temp;
				}
			}

	for (i = 0; i < pdim; i++)
		for (j = 0; j < pdim; j++)
			S[i * pdim + j] = S[i * pdim + j];

    /*
	printf("the eigenvalue is :\n");
	for (i = 0; i < pdim; i++)
		printf("%f ", matrix[i * (i + 1) / 2 + i]);
	printf("\n");

	printf("the S matrix:\n");
	for (i=0; i<pdim; i++){
		for (j=0; j<pdim; j++)
			printf("%f, ",S[i*pdim+j]);
		printf("\n");
	}	
    */


	/* orientation sign (old) */
	for (i=0; i<pdim; i++){
		if (S[i] < 0)
			for (j=0; j<pdim; j++)
				S[j*pdim+i] = -S[j*pdim+i];		
	}

    /* orientation sign (new) 
    //get the maximum variation dimension of the last eigenvector (the one with largest eigen value)
    double maxvar=0, tempf;
	int   maxi=-1;
	for (i=0; i<pdim; i++)
	{
        tempf=fabs(S[pdim-1+i*pdim]);
	    if (fabs(S[pdim-1+i*pdim])>maxvar)
		{
           maxvar=fabs(S[pdim-1+i*pdim]);
		   maxi=i;
		}
	}
    printf("maxvar:%f, maxi:%d\n", maxvar, maxi);
	for (i=0; i<pdim; i++){
		if (maxi<3 && S[i+2*pdim] < 0) 
			for (j=0; j<pdim; j++)
				S[j*pdim+i] = -S[j*pdim+i];	
		else if (maxi>=3 && maxi<6 && S[i+5*pdim] < 0)
			for (j=0; j<pdim; j++)
				S[j*pdim+i] = -S[j*pdim+i];	
		else if (maxi>=6 && S[i+8*pdim] < 0)
			for (j=0; j<pdim; j++)
				S[j*pdim+i] = -S[j*pdim+i];
	}*/

	for (i=0; i<pdim; i++)
		for (j=0; j<pdim; j++)
			fprintf(pofp, "%f ",S[j*pdim+i] * 255);
	fprintf(pofp, "\n");	
}



/*******************************EM Clustering**************************************/
void EM::EMCluster(char *pinfile, char *poutfile, char *porderfile, float pepsilon, float pdistepsilon, float pmindensity){
	
	/* file pointer to the data file */
	FILE *ifp;
	FILE *ofp;

	/* the maximum liklyhood */
	float epsilon;

	/* iteration variable */
	int i, j, k, m;

	/* the number of points process as so far */
	float NUM;

	/* stopping criterion for RS set */
	float E, newE;

	/* local variable */
	float tempresult;
	float *tempPRlx;
		
	float *tempv;
	float PRx;
	int count;
	int sigma_index;
		
	tempv = (float *) malloc (sizeof (float) * dim);
	PRxl = (float *) malloc (sizeof (float) * numcluster);
	PRlx = (float *) malloc (sizeof (float) * samplesize * pagecap);
	uXXT = (float *) malloc (sizeof (float) * numcluster * dim * (dim+1)/2);
	uX = (float *) malloc (sizeof (float) * numcluster * dim);
	dist = (float *) malloc (sizeof (float) * numcluster * (numcluster-1)/2);
	tempPRlx = (float *) malloc (sizeof (float) * numcluster);

	/* memory allocation test */
	if (tempv == NULL){
		printf("Not enough memory for tempv!!\n");
		exit(0);
	}
	if (PRxl == NULL){
		printf("Not enough memory for PRxl!!\n");
		exit(0);
	}
	if (PRlx == NULL){
		printf("Not enough memory for PRlx!!\n");
		exit(0);
	}
	if (uXXT == NULL){
		printf("Not enough memory for uXXT!!\n");
		exit(0);
	}
	if (uX == NULL){
		printf("Not enough memory for uX!!\n");
		exit(0);
	}
	if (dist == NULL){
		printf("Not enough memory for dist!!\n");
		exit(0);
	}
	if (tempPRlx == NULL){
		printf("Not enough memory for tempPRlx!!\n");
		exit(0);
	}
		
	/* initialization */
	infile = pinfile;
	outfile = poutfile;
	epsilon = pepsilon;
	
	E = 0;
	newE = 100000;
	NUM = 0;

	
	int tmp_HowMany = samplesize;
	int t, tmp_i, tmp_j, tmp_k;
	float tmp_f;
	tmp_ItemList = (struct S_BigEnoughItem*)malloc( pagecap * sizeof(struct S_BigEnoughItem) );
	if (tmp_ItemList == NULL){
		printf("Not enough memory for tmp_ItemList!\n");
		exit(0);
	}

		
	/* open the file for recording cluster information */
	if ((ofp = fopen(poutfile, "w")) == NULL)
	{
		printf("Outfile not opened!");
		exit(0);
	}

	//initialization for the sumprob procedure, EM preparation

//	Timer emtime;
//	emtime.start();
	int iterationnum=0;

	/***********************the Outer Loop******************************************/
	//if the stop criterion for outer loop is not satisfied
	float g_Factor = float(pow(2 * 3.14, dim));
	if (g_Factor > 1000000)
		g_Factor = 10000000;

	/* the outer iteration loop, begins */
	while ((fabs(E - newE) > epsilon) && (iterationnum < iterationtime)){
		iterationnum = iterationnum + 1;
		printf("%dth iteration begin...\n", iterationnum);
		E = newE;
		newE = 0;
		
		/* initialization for inner loop */
		for (i = 0; i < numcluster; i++){
			N[i] = 0;			  
			for (j = 0; j < dim; j++){
				uX[i * dim + j] = 0;
				for (k = 0; k <= j; k++)
					uXXT[i * dim * (dim+1)/2 + j * (j+1)/2 + k] = 0;				  
			}
		}
		
		//End of M-Step;???
		//determinant of sigma;
		for (j = 0; j < numcluster; j++){
			determinant[j] = Det(&sigma[j * sigmasize], dim);
//			if (iterationnum != 1)
//				printf("the determinant of cluster %d is %f;\n", j, determinant[j]);
			if (iterationnum == 1)
				determinant[j] = float(sqrt(determinant[j]));
			else
				determinant[j] = float(sqrt(g_Factor * fabs(determinant[j])));
//			if (iterationnum != 1)
//				printf("the determinant of cluster %d after modify is %f;\n", j, determinant[j]);
		}
					 
		/* inverse of sigma */
		InverseMatrix(sigma, dim, numcluster);

		count = 0;

		/* open the data file */	
		if ((ifp = fopen(infile, "r")) == NULL)
		{
			printf("Infile Not Opened!");		
			exit(0);
		}

		/* the first iteration */
		if (iterationnum == 1)
			while (!feof(ifp)){	

				/* fscan each data point in the file */ 
				for (j = 0; j < dim; j++){
					fscanf(ifp, "%f\t", &RS[j]);
				}
				fscanf(ifp, "\n");
						
        /***********************E-step for Singleton Records***************************/
		
			    /* calculate the Pr(x|l) for each cluster and store it in PRxl[l] */
				for (j = 0; j < numcluster; j++){
					PRxl[j] = 0;
					/* x-ul (ul: the vector of cluster l) */
					for (k = 0; k < dim; k++)
						tempv[k] = RS[k]-u[j * dim + k];
				
					/* (x-ul)*inversesigma*(x-ul)T (the mahalanobis distance) */
					/* tempv: 1*dim; inversesigma: dim(m)*dim(k) */
					for (k = 0; k < dim; k++){
						tempresult = 0;
						for (m = 0; m < dim; m++){
							if (k < m)
								sigma_index = j * sigmasize + m * (m + 1) / 2 + k;
							else 
								sigma_index = j * sigmasize + k * (k + 1) / 2 + m;
	
								tempresult = tempresult + tempv[m] * sigma[sigma_index];
						}
					
						PRxl[j] = PRxl[j] + tempresult * tempv[k];
					}
				
					PRxl[j] = float(fabs(PRxl[j]));    
				
				}

				/* calculate PRx */
				PRx = 0;
				for (j = 0; j < numcluster; j++){
					if (PRxl[j] < 600)
						PRxl[j] = float(exp(- PRxl[j] / 2)/ determinant[j]);
					else
						PRxl[j] = 0;
					}
				for (j = 0; j < numcluster; j++){
					PRx = PRx + PRxl[j];
				}				

				/* update newE */
				if (PRx > 0) newE = newE + float(log(PRx));
			
				for (j = 0; j < numcluster; j++){
					/* calculate PRlx */
					if (PRx > 0) tempPRlx[j] = PRxl[j] / PRx;
					else tempPRlx[j] = 0;

					/* uX */
					for (k = 0; k < dim; k++){
						tempv[k] = RS[k]-u[j * dim + k];
						uX[j * dim + k] = uX[j * dim + k] + RS[k] * tempPRlx[j];
					}
				
					N[j] = N[j] + tempPRlx[j];

				    /* uXXT */
					for (k = 0; k < dim; k++){
						for (m = 0; m <=k; m++)
							uXXT[j * dim * (dim+1)/2 + k * (k+1)/2 + m] = uXXT[j * dim * (dim+1)/2 + k * (k+1)/2 + m] +
							   tempv[k] * tempv[m] * tempPRlx[j];					
					}
				}
			
				/* set up an itemlist for each data point */
				/* store the nearest 'samplesize' cluster ids */
				tmp_ItemList[count].ItemNO = samplesize;//samplesize: compressed number of cluster
				tmp_ItemList[count].ItemList = (short *)malloc( samplesize * sizeof(short) );
				if (tmp_ItemList[count].ItemList == NULL){
					printf("Not enough memory!!\n");
					exit(0);
				}

				tmp_i = 0;
				for ( j = 0; j < numcluster; j++ ){
					tmp_f = tempPRlx[j];

					/* if the itemlist is not full */
					if ( tmp_i < tmp_HowMany ){
						tmp_i++;
						tmp_ItemList[count].ItemList[tmp_i - 1] = -1;			
					}	
					
					/* compare the distance of the data to the current cluster to that of other clusters */
					/* if it is nearer than ONE of them, then insert it and move it */
					for ( k = tmp_i - 1; k >= 0; k-- ){
						if ( tmp_ItemList[count].ItemList[k] != -1 && 
							//tmp_f <= tempPRlx[tmp_ItemList[count].ItemList[k] ] ) 
							tmp_f > tempPRlx[tmp_ItemList[count].ItemList[k] ] ) 
							break;
					}
					/////////////////////////////////////////////tmp_f ERROR

					if ( k < tmp_HowMany-1 ){
						for ( t = tmp_i - 2; t > k; t-- ){
							tmp_ItemList[count].ItemList[t + 1] = tmp_ItemList[count].ItemList[t];
//							PRlx[count *  samplesize + t+1] = PRlx[count * samplesize + t];
						}
						tmp_ItemList[count].ItemList[k + 1] = j;
					}
				}

				for (j = 0; j < samplesize; j++){
					tmp_i = tmp_ItemList[count].ItemList[j];
//					printf(" tmp_ItemList[%d].ItemList[%d] = %d\n", count, j,  tmp_ItemList[count].ItemList[j]);
					PRlx[count*samplesize+j] = tempPRlx[tmp_i];
				}
            
				/* count: the index of the current processed data point */
				count = count + 1;
			}//while !eof ifp??
		else { //when it is not the first iteration

            /* for each point in the file */
			while (!feof(ifp)){	
				 
				for (j = 0; j < dim; j++){
					fscanf(ifp, "%f\t", &RS[j]);
				}
				fscanf(ifp, "\n");
						
        /***********************E-step for Singleton Records***************************/
		        
				/* calculate the Pr(x|l) for the samplesize clusters and store it in PRxl[l] */
				int tmp_cluster;
				for (j = 0; j < samplesize; j++){
					tmp_cluster = tmp_ItemList[count].ItemList[j];
					PRxl[j] = 0;

					/* x-ul */
					for (k = 0; k < dim; k++){						
						tempv[k] = RS[k]-u[tmp_cluster * dim + k];
					}
				
				
					/* (x-ul)*inversesigma*(x-ul)T (the mahalanobis distance) */
					for (k = 0; k < dim; k++){
						tempresult = 0;
						for (m = 0; m < dim; m++){
							if (k < m)
								sigma_index = tmp_cluster * sigmasize + m * (m + 1) / 2 + k;
							else 
								sigma_index = tmp_cluster * sigmasize + k * (k + 1) / 2 + m;
	
								tempresult = tempresult + tempv[m] * sigma[sigma_index];
						}
					
						PRxl[j] = PRxl[j] + tempresult * tempv[k];
					}

					PRxl[j] = float(fabs(PRxl[j]));  
					if (PRxl[j] < 600)
						PRxl[j] = float(exp(- PRxl[j] / 2) / determinant[tmp_cluster]);
					else
						PRxl[j] = 0; 
				}			

				/* calculate PRx */	
				PRx = 0;	
				for (j = 0; j < samplesize; j++){
					PRx = PRx + PRxl[j];
				}				

				/* update newE */
				if (PRx > 0) newE = newE + float(log(PRx)/100);
			
				/* calculate PRlx */
				for (j = 0; j < samplesize; j++){
					if (PRx > 0) PRlx[count * samplesize+j] = PRxl[j] / PRx;
					else PRlx[count*samplesize+j] = 0;

					tmp_cluster = tmp_ItemList[count].ItemList[j];

					/* uX */
					for (k = 0; k < dim; k++){
						tempv[k] = RS[k]-u[tmp_cluster * dim + k];
						uX[tmp_cluster * dim + k] = uX[tmp_cluster * dim + k] + RS[k] * PRlx[count*samplesize+j];
					}
				
					N[tmp_cluster] = N[tmp_cluster] + PRlx[count*samplesize+j];

				    /* uXXT */
					for (k = 0; k < dim; k++){
						for (m = 0; m <=k; m++)
							uXXT[tmp_cluster * dim * (dim+1)/2 + k * (k+1)/2 + m] = uXXT[tmp_cluster * dim * (dim+1)/2 + k * (k+1)/2 + m] +
							   tempv[k] * tempv[m] * PRlx[count * samplesize + j];					
					}
				}
			
				count = count + 1;
			}		
		}//end of else
		

		fclose(ifp);

	/************************************* End of E-Step;******************************/

    /************************************* M step begins ******************************/


		for (j = 0; j < numcluster; j++){
			if (N[j] + 1 == 1)
				InverseMatrix(&sigma[j * sigmasize], dim, 1);
			else {
				for (k = 0; k < dim; k++){
					/* update u */
					u[j * dim + k] = uX[j * dim + k] / N[j];

					/* update sigma */
					for (m = 0; m < dim; m++){
						if (k < m)
							sigma_index = j * sigmasize + m * (m + 1) / 2 + k;
						else 
							sigma_index = j * sigmasize + k * (k + 1) / 2 + m;

						sigma[sigma_index] = uXXT[sigma_index] / N[j];
						if (fabs(sigma[sigma_index]) < 0.001)
							if (k == m) sigma[sigma_index] = float(0.001);						
							else sigma[sigma_index] = 0;						
					}
				}
			}
		}
							
	}//end of the iteration outloop and fabs(E-newE)<

    /* newly added */
    /* calculate the average number of points that belong to one microcluster */
    float sumprobability=0;
	//float temp;
	for ( j = 0; j < samplesize; j++ )
        for ( i = 0; i < pagecap; i++ )		
		     sumprobability=sumprobability+PRlx[ i * samplesize + j ];
    printf("\nAverage Data Points for a microcluster: %.2f\n\n", sumprobability/samplesize);

	//newly added
	/* record point probability */
	FILE *tt;
	int member;
	double maxmembership;
	tt=fopen("membership.txt", "w");
	for (i = 0; i < pagecap; i++)
	{
		member=-1;
		maxmembership=0;
		for (j = 0; j < samplesize; j++)
			if (PRlx[i*samplesize+j]>maxmembership)
			{
				maxmembership=PRlx[i*samplesize+j];
				member=tmp_ItemList[i].ItemList[j];
			}
			//fprintf(tt, "%d (%f), ", tmp_ItemList[i].ItemList[j], PRlx[i*samplesize+j]);
		fprintf(tt, "%d\n", member);
    }
	fclose(tt);

	/* initialize the distance array */
	for (i = 0; i < numcluster; i++)
		for (j = 0; j < i; j++)
			dist[i * (i-1)/2 + j] = 0;

    printf( "distance computation begins ... \n" );
    //emtime.start();
	int dist_index, tmp_SampleSize;	

	tmp_SampleSize = samplesize;
	tmp_HowMany = samplesize * compressratio2;
	int* tmp_Array;
	tmp_Array = (int *) malloc (sizeof(int) * samplesize);

	/* only keep pagecap data points  */
	for ( i = 0; i < pagecap; i++ )
	{
		tmp_i = 0;
		for ( j = 0; j < tmp_SampleSize; j++ )
		{
			tmp_f = PRlx[ i * tmp_SampleSize + j ];
			if ( tmp_i < tmp_HowMany )
			{
				tmp_i++;
				tmp_Array[tmp_i - 1] = -1;
			}

			for ( k = tmp_i - 1; k >= 0; k-- )
			{
				//if ( tmp_Array[k] != -1 && tmp_f <= PRlx[ i * tmp_SampleSize + k ] )
				if ( tmp_Array[k] != -1 && tmp_f > PRlx[ i * tmp_SampleSize + k ] )

					 break;
			}
		///////////////////////////////ERROR tmp_f

			if ( k < tmp_HowMany - 1 )
			{
				for ( t = tmp_i - 2; t > k; t-- )
				{
					tmp_Array[t + 1] = tmp_Array[t];
					PRlx[ i * tmp_SampleSize + t + 1] = PRlx[ i * tmp_SampleSize + t ];
				}
				tmp_Array[k + 1] = tmp_ItemList[i].ItemList[j];
				PRlx[ i * tmp_SampleSize + k + 1] = tmp_f;
			}
		}

		// update the top tmp_HowMany elements of tmp_Itemlist[i]
		for ( k = 0; k < tmp_HowMany; k++ )
		{
			tmp_ItemList[i].ItemList[k] = tmp_Array[k];			
		}

	}//end of i<pagecap
   

	free(tmp_Array);

	/* calculate those possiblities above threshold (first)*/
	for (i = 0; i < pagecap; i++)
	{	
		for (j = 0; j < tmp_HowMany; j++)
		{
			tmp_j = tmp_ItemList[i].ItemList[j];
			
			for (k = j + 1; k < tmp_HowMany; k++)
			{				
				tmp_k = tmp_ItemList[i].ItemList[k];
				if (tmp_k > tmp_j)
					dist_index = tmp_k * (tmp_k-1)/2 + tmp_j;
				else 
					dist_index = tmp_j * (tmp_j-1)/2 + tmp_k;

				// update the distance matrix
				dist[dist_index] = dist[dist_index] + PRlx[i*samplesize+j] * PRlx[i*samplesize+k];				
			}
		}
	} 


    /* show distance matrix 
	float nearestdist;
	int   neighbor;
	printf("distance matrix....\n");
	for (i = 0; i < numcluster; i++)// i=0; i < numcluster
	{
		nearestdist=0;
		for (j = 0; j < i; j++)			
		{
			if (dist[i * (i-1)/2 + j]>nearestdist)
			{
                nearestdist=dist[i * (i-1)/2 + j];
				neighbor=j;
			}
		}
		if (nearestdist>0)
		    printf("dist(%d, %d): %f\t", i, neighbor, nearestdist);

        nearestdist=0;
		for (j = i+1; j < numcluster; j++)			
		{
			//printf("dist(%d, %d): %f ",i, j, dist[i * (i-1)/2 + j]);
			//if (j%3==2) printf("\n");
			if (dist[j * (j-1)/2 + i]>nearestdist)
			{
                nearestdist=dist[j * (j-1)/2 + i];
				neighbor=j;
			}
		}
		if (nearestdist>0)
		   printf("dist(%d, %d): %f\n", i, neighbor, nearestdist);
	}*/
     
  
	for ( i = 0; i < pagecap; i++ ) 
		free( tmp_ItemList[i].ItemList );
	free( tmp_ItemList );
	free(PRlx);	
//emtime.end();
//printf( "compute distance end\n" );

	/**********************************************************************************/
	/**********************************end of EM algorithm*****************************/
	/**********************************************************************************/

//	emtime.start();
	optics(porderfile, pdistepsilon, pmindensity);
//	emtime.end();


	free(sigma);
 	//free(PRlx);
	free(uXXT);
	free(uX);
	free(determinant);
	free(RS);
	free(tempv);
	free(PRxl);
 	free(dist);
	free(tempPRlx);


	//open the ordered file for read
	if ((ifp = fopen(infile, "r")) == NULL)
		printf("File Not Opened!");
		  
	float tempfval;
	//purge the data remained in the RSset
	for (i = 0; i < pagecap; i++){
		for (j = 0; j < dim; j++){
			fscanf(ifp,"%f ", &tempfval);
			fprintf(ofp, "%f ", tempfval);
		}
		fscanf(ifp, "\n");
		fprintf(ofp, "\n");
	}

	fclose(ifp);		
	
	if ((ifp = fopen(porderfile,"r")) == NULL)
	{
		printf("File Not Opened!");
		exit(0);
	}

	float reach;
	int index,temp1,temp2;

	for (i = 0; i < numcluster; i++){
		//move to the last eigenvector
		fscanf(ifp, "%f %d ", &reach, &index);
		for(temp1=1;temp1<dim;temp1++)
			for(temp2=0;temp2<dim;temp2++)
               fscanf(ifp, "%f ", &reach);
        fscanf(ifp, "\n");

		for (j = 0; j < dim; j++)
			fprintf(ofp, "%f ", u[index * dim +j]);
		for (j = 0; j < dim; j++){
			fscanf(ifp, "%f ", &reach);
			fprintf(ofp, "%f ", reach);
		}
		//for (j = 0; j < (dim * dim-dim); j++)
			//fscanf(ifp, "%f ", &reach);
		
		fprintf(ofp, "%d\n", index);
		
	}

	free(u);

	fclose(ifp);
	fclose(ofp);

}

void EM::filewrite(char *poutfile, float *pu, float *psigma, int pdim, int pnumcluster)
{
	FILE *ofp;
	int i;
	
	//format: u vector of the cluster;
	//        sigma matrix of the cluster;

	if ((ofp = fopen(poutfile, "wb")) == NULL)
	{
		printf("File Not Opened!");
		exit(0);
	}

	for (i = 0; i < pnumcluster; i++){
		fwrite(&pu[i * pdim], sizeof(float), pdim, ofp);
		fwrite(&psigma[i * pdim * pdim], sizeof(float), pdim * pdim, ofp);
	}

	fclose(ofp);
}

float EM::probdifference(int pnumcluster1, int pnumcluster2)
{
	double probdif;
	int i;

	probdif=0;
	for (i = 0; i < pagecap; i++){
		probdif = probdif + fabs(PRlx[i * numcluster + pnumcluster1] * PRlx[i * numcluster + pnumcluster2])*100;
	}

	return float(probdif);
}

/* a subroutine called by optic */
/* pclusternum: current microcluster */
/* pprocessed: indicator array */
void EM::expandclusterorder(int pclusternum, float pdisteps, float pmindensity, int *pprocessed, FILE *pofp)
{
	int *neighindex;
	float *neighdif;
	int j;
	int neighcount;
	int neighcluster;
	int index1, index2;
	
	neighindex = (int *) determinant;
	neighdif = PRxl;

	/* update the processing record */
	pprocessed[pclusternum] = 1;
	processedcluster = processedcluster + 1;
//	printf("%d clusters processed\n", processedcluster);
	cout << "\r\t"<<(int)(processedcluster*100.0/numcluster) << "% [" 
			<< processedcluster << " of " << numcluster << " microclusters]" << flush;

	neighcount = neighbour(pclusternum, pdisteps, neighindex, neighdif, pprocessed);

	fprintf(pofp, "%f %d ", -1.000, pclusternum);
//	printf("cluster %d:\n", pclusternum);

	/* record the eigen value of current microcluster to file */
	eigenvalue(pofp, &sigma[pclusternum * sigmasize], dim);

	/* choose the nearest neighbour as the next move, then choose its nearest neighbour ... */
	while (neighcount > 0){
		neighcluster = neighindex[0];
		pprocessed[neighcluster] = 1;
//		printf("the neighbor of cluster %d is %d\n", pclusternum, neighcluster);
//		for (j = 0; j < pagecap; j++)
//			PRlx[j * numcluster + pclusternum] = (PRlx[j * numcluster + j] + PRlx[j * numcluster + neighcluster]);
		
		/* update the distance matrix */
		for (j = 0; j < numcluster; j++)
			if ((pprocessed[j] ==0) && (j!=pclusternum)){
				if (j < pclusternum){
					index1 = pclusternum * (pclusternum-1)/2 + j;
				}
				else{
					index1 = j*(j-1)/2 + pclusternum;
				}
				if (j < neighcluster){
					index2 = neighcluster * (neighcluster-1)/2 + j;
				}
				else{
					index2 = j*(j-1)/2 + neighcluster;
				}
				/*sum
				dist[index1] = dist[index2] + dist[index1];*/

				/*maximum*/
				if (dist[index2]>dist[index1])
					dist[index1]=dist[index2];
				
			}

//			PRlx[pclusternum * pagecap + j] = (PRlx[pclusternum * pagecap + j] + PRlx[neighcluster * pagecap + j]);
//		for (j = 0; j < pagecap; j++){
//			PRlx[j * numcluster + pclusternum] = 4*pow(PRlx[j * numcluster + pclusternum],2) + 4*pow(PRlx[j * numcluster +neighcluster],2) - pow(PRlx[j* numcluster + neighcluster]*PRlx[j*numcluster+pclusternum],2);
//			if (PRlx[j*numcluster+pclusternum] < 0){
//				printf("the new prob is of point %d: %f\n", j, PRlx[j*numcluster+pclusternum]);
//				PRlx[j*numcluster+pclusternum] = 0;
//			}
//			else 
//				PRlx[j*numcluster+pclusternum] = sqrt(PRlx[j*numcluster+pclusternum])/2;
//		}
//		for (j = 0; j < pagecap; j++)
//			if (PRlx[j * numcluster + pclusternum] < PRlx[j * numcluster +neighcluster])
//				PRlx[j * numcluster + pclusternum] = PRlx[j * numcluster +neighcluster];

		fprintf(pofp, "%f %d ", neighdif[0], neighcluster);
//		printf("cluster %d:\n", neighcluster);
		eigenvalue(pofp, &sigma[neighcluster * sigmasize], dim);
		processedcluster = processedcluster + 1;
//		pclusternum = neighcluster;
		cout << "\r\t"<<(int)(processedcluster*100.0/numcluster) << "% [" 
			<< processedcluster << " of " << numcluster << " microclusters]" << flush;
//		printf("cluster %d processed\n",neighcluster);
//		printf("%d clusters processed\n", processedcluster);
		
		/* reset neighindex and neighdif */
		for (j = 0; j < numcluster; j++){
			neighindex[j] = -1;
			neighdif[j] = -1;
		}
		neighcount = neighbour(pclusternum, pdisteps, neighindex, neighdif, pprocessed);
		//neighcount = neighbour(neighcluster, pdisteps, neighindex, neighdif, pprocessed);
			
	}
}


void EM::optics(char *poutfile, float pdistepsi, float pmindensity)
{
	//flag to indicate whether this cluster is processed or not
	int *processed;
	FILE *ofp;                      //file for output
	int i;
	

	/* processed: indicator array, indicating whether the point is processed or not */
	processed = (int *) N;
	processedcluster = 0;
	
	for (i = 0; i < numcluster; i++){
		processed[i] = 0;		
	}

	/* open ordered file for write */
	if ((ofp = fopen(poutfile,"w")) == NULL)
	{
		printf("File Not Opened!");
		exit(0);
	}

	/* check each microcluster */
	for (i = 0; i < numcluster; i++){
		//if not processed
		if (processed[i] == 0)
				expandclusterorder(i, pdistepsi, pmindensity, processed, ofp);
/*			else {
				fprintf(ofp, "%f %d ", -1.000, i);
//				printf("cluster %d:\n", i);
				eigenvalue(ofp, &sigma[i * sigmasize], dim);
				processed[i] = 1;
				processedcluster = processedcluster + 1;
				cout << "\r\t"<<(int)(processedcluster*100.0/numcluster) << "% [" 
					<< processedcluster << " of " << numcluster << " microclusters]" << flush;
//				printf("cluster %d processed\n",i);
//				printf("%d clusters processed\n", processedcluster);
			}*/
	}
		
	//close ordered file
	fclose(ofp);
}

/* find the neighbours of the current microcluster "clusternum" */
int EM::neighbour(int clusternum, float pdistepsi, int *pneighindex, float *pneighdif, int *pprocessed)
{
	int i,j;
	float truedif, tempdif;
	int neighcount;
	int tempindex;

	neighcount = 0;

	/* record the microclusters within distance of pdistepsi */
	for (i = 0; i < numcluster; i++)
		if ((pprocessed[i] == 0)){
			if (i < clusternum)
				truedif = dist[clusternum*(clusternum-1)/2+i];
			else
				truedif = dist[i*(i-1)/2+clusternum];
			if ((truedif > pdistepsi)){
						pneighindex[neighcount] = i;
						pneighdif[neighcount] = truedif;
						neighcount = neighcount + 1;
					}
		}


	/* order the neighbor list */
	for (i = 0; i < neighcount; i++)
		for(j = 0; j < neighcount -1 - i; j++)
			if (pneighdif[j] < pneighdif[j + 1]){
				tempdif = pneighdif[j];
				tempindex = pneighindex[j];
				pneighdif[j] = pneighdif[j + 1];
				pneighindex[j] = pneighindex[j + 1];
				pneighdif[j + 1] = tempdif;
				pneighindex[j + 1] = tempindex;
			}

//	for (i = 0; i < neighcount; i++)
//		printf("the %dth neighbor of cluster %d is cluster %d: %f\n", i, clusternum, pneighindex[i], pneighdif[i]);
				
//           	printf("neighbor searching end...\n");
	return(neighcount);

}

int EM::updatedif(int *pseedflag, int *pcluster2seed, int *pseedlist, float *pseeddif, int *pneighindex, float *pneighdif, int pseedcount, int pneighcount)
{
	int i;
	int seedindex;
	float newreach;

	//update seedlist
	for (i = 0; i < pneighcount; i++)
		//the cluster is not in the seedlist
		if (pseedflag[pneighindex[i]] == 0){
			pseedlist[pseedcount] = pneighindex[i];
			pseeddif[pseedcount] = pneighdif[i];
			pseedflag[pneighindex[i]] = 1;
			pcluster2seed[pneighindex[i]] = pseedcount;
			pseedcount = pseedcount + 1;			
		}
		//the cluster is already in the seedlist;
		else {
//			if (pneighdif[i] > pcoredist)
				newreach = pneighdif[i];
//			else 
//				newreach = pcoredist;

			seedindex = pcluster2seed[pneighindex[i]];
			if (pseeddif[seedindex] < newreach)
				pseeddif[seedindex] = newreach;
		}

	//reorder the seedlist
	orderseedlist(pseedcount, pseedlist, pseeddif, pcluster2seed);

	return pseedcount;

}

void EM::orderseedlist(int pseedcount, int *pseedlist, float *pseeddif, int *pcluster2seed)
{
	int i, j;
	int tempindex, oldindex;
	float olddif;

	for(i = 0; i < pseedcount; i++){
		for (j = 0; j < pseedcount - i -1; j++)			
			if (pseeddif[j] < pseeddif[j + 1]){
				oldindex = pseedlist[j];
				olddif = pseeddif[j];
				pcluster2seed[pseedlist[j + 1]] = j;
				pseedlist[j] = pseedlist[j + 1];
				pseeddif[j] = pseeddif[j + 1];
				tempindex = j + 1;
				pseedlist[j + 1] = oldindex;
				pseeddif[j + 1] = olddif;
				pcluster2seed[oldindex] = j + 1;
			}		
	}
}


void EM::popseedlist(int *pseedlist, float *pseeddif, int *pcluster2seed, int seedcount)
{
	int k;
	for (k = 1; k < seedcount; k++){
		pseedlist[k-1] = pseedlist[k];
		pseeddif[k-1] = pseeddif[k];
		pcluster2seed[pseedlist[k]] = pcluster2seed[pseedlist[k]] - 1;
	}
}

float EM::setcoredistance(float pN, float pmindensity, int pneighcount, float *pneighdif, int *pneighindex)
{
	printf("setting core distance...\n");
	float* distorder;
	float coredistance;
	float temp;
	int i, j;;
	float density;
	int *distindex;
	int tempindex;

	density = pN;
	for (i = 0; i < pneighcount; i++)
		density = density + N[pneighindex[i]];


	if (density < pmindensity)
		/* the sum of density does not satisfy the minimum density requirement */
		return -1;
	else{
		distorder = (float *) malloc (sizeof(float) * pneighcount);
		distindex = (int *) malloc (sizeof(int) * pneighcount);

		//initial the distorder array with the pminpts values at the begining of the pneighdif array;
		for (i = 0; i < pneighcount; i++){
			distorder[i] = pneighdif[i];
			distindex[i] = pneighindex[i];
		}

		/* order the initial distorder array in ascending order */
		for (i = 0; i < pneighcount; i++){
			for (j = 0; j < pneighcount - i - 1; j++)
				if (distorder[j] > distorder[j + 1]){
					temp = distorder[j];
					tempindex = distindex[j];
					distorder[j] = distorder[j + 1];
					distindex[j] = distindex[j + 1];
					distorder[j + 1] = temp;		
					distindex[j + 1] = tempindex;
				}				
		}

		/* locate the first time when the minimum density is satisfied */
		density = pN;
		for (j = 0; j < pneighcount; j++){
			density = density + N[distindex[j]];
			if (density >= pmindensity) break;
		}
		coredistance = distorder[j];
		free(distorder);
		free(distindex);
		return coredistance;
	}

	printf("core distance setting finished\n");

}

float EM::Det(float *matrix, int pdim){
	double determinant;
	float *tempmatrix;
	int flag;
	int i, j, k;
	int det_index;

	flag = 1;
	
	tempmatrix = tempmem;

	/* assign values to tempmatrix */
	for (i = 0; i < pdim; i++)
		for (j = 0; j < pdim; j++){
			if (i < j)
				det_index = (j + 1) * j / 2 + i;
			else 
				det_index = (i + 1) * i / 2 + j;

			tempmatrix[i * pdim + j] = matrix[det_index];
		}

	/* LU decomposition */
	for (i = 0; i < pdim-1; i++){
		if (fabs(tempmatrix[i * pdim + i])+1.0 == 1.0){
			flag = 0;
			determinant = 0.00001;
			return float(determinant);
		}
		/* divide the elements of column i row >i by the element at [i,i] */
		for (j = i + 1; j< pdim; j++)
			tempmatrix[j * pdim + i] = tempmatrix[j * pdim + i] / tempmatrix[i * pdim + i];
		for (j = i + 1; j < pdim; j++)
			for (k = i + 1; k < pdim; k++)
				tempmatrix[j * pdim + k] = tempmatrix[j * pdim + k] - tempmatrix[j * pdim + i] * tempmatrix[i * pdim + k];
	}

	determinant = 1;
	for (i = 0; i < pdim; i++)
		determinant = determinant * tempmatrix [i * pdim + i];

	if (determinant + 1.0 == 1.0)
		determinant = 0.000001;
	return float(determinant);	
}