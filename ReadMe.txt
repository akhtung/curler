****************** Parameter Setting of Curler Algorithm ******************
datafile:           the dataset file, with no class labels, one line for one data.
#cluster:           the number of microclusters.
#dimension:         the dimensionality of the dataset.
#data:              the number of data points.
epsilon_likelihood: probability log likelihood threshold.
epsilon_coshare:    neighborhood co-sharing level threshold (default 0).
microclus file:     the microcluser file of #cluster randomely sampled data
MaxLoopNum:         maximum number of iteration times (2-15 is sufficient). 
basedim:            basedim value of orientation vectors is kept positive (default 0, range [0,#dimension-1]).
       


******************     Output of Curler Algorithm     ******************
output_em.txt: first part:      the original data 
                               (the d dimension values of each data)
               second part:     the microclusters.
                               (the d dimension values of mean point )
                               (one eigenvector of d dimensions)
                               (the index of the microcluster)
 
output_expand.txt:the ordered microclusters, one line for one microcluster
                  ( NND nearest neighbour distance) 
                  ( nearest microcluster neighbour index (may be not consecutive))
                  ( d eigen vectors, each of which has d dimension values )
 	          
membership.txt:  the index of the microcluster to which each data belongs to.



******************     Running Example     ******************
iris dataset:
Step 1 (dos): clustering
Curler iris.txt 150 4 150 0 0 iris_seed150.txt 8 0
Step 2 (matlab): visualization
NNCO(150,4,'output_expand.txt');

image dataset:
Step 1 (dos): clustering
Curler image.txt 500 16 2310 0 0 image_seed500.txt 8 0
Step 2 (matlab): visualization
NNCO(500,16,'output_expand.txt');

cubic dataset:
Step 1 (dos): clustering
Curler cubic.txt 150 2 589 0 0 cubic_seed150.txt 3 0
Step 2 (matlab): visualization
NNCO(150,2,'output_expand.txt');