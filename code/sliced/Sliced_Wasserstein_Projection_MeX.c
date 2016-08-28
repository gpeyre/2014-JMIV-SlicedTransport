#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "math.h"
#include "mex.h"

/*
 * (function call in Matlab : [X,E] = Sliced_Wasserstein_Projection_MeX(X0, Y, step_descent, Dir, niter,opt_stochastic, opt_hessian, opt_verbose);) 
 *
 * Inputs : X0 = input point clouds [d x N] matrix of N vectors of dimension d
 *          Y  = target opint clouds
 *          step_descent = descent parameter (should be <1)
 *          Dir = either the number of directions or a [d x k] matrix of k vectors
 *          niter = maximum number of iterations
 *          opt_stochastic = option to use random directions at each step (stochastic gradient descent)
 *          opt_hessian = option to use hessian normalization (Newton descent)
 *          opt_verbose = option to verbose
 * 
 * Outputs : X = ouput point-cloud (X-X0 is the sliced wasserstein transport flow)
 *           E = energy monitoring
 *
 * Exemple 1 : 
              X0 = randn(2,1e5); Y = [2,0;0,1/2]*X0+2;
              tic, [X,E] = Sliced_Wasserstein_Projection_MeX(X0,Y,1,1e1,1e2,0,0,0), toc
              figure, plot(X0(1,:),X0(2,:),'ob'), hold on, 
              plot(Y(1,:),Y(2,:),'or'), plot(X(1,:),X(2,:),'xg')
 *
 * Julien.Rabin@unicaen.fr (c) 2013 - GREYC, University of Caen
 *
*/


#define DEBUG 0
#define USE_PREVIOUS_ORDERING 1

typedef struct {
    int D; // dimension
    int N; // number of points
    double step; // step parameter for gradient descent
    int ndir; // number of directions (slices)
    int niter; // number of iterations
    double * Dir; // array of directions ( D x ndir)
    int flag_stochastic; // random direction at each step
    int flag_hessian; // Hessian normalization (Newton descent)
    int verbose;
} OPT;


// prototypes of auxillary functions
mxArray * Sliced_Wasserstein_Projection(mxArray *X,const mxArray *Y,OPT options);
double rand_normal(double mean, double stddev);
void copy_array(double *in,double *out,int N);
// permute from index array cast as 'double'
void self_permute_array(double * Z, double * z, double * p, int nrow, int ncol);
void perform_permute_array(double * in, double * out, double * p, int ndir, int N) ;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/* - - - - - - - -
 *  Declarations
 * - - - - - - - -
 */ 

    mxArray *X0Data, *YData;
    //double * X0, *X, *Y;
    double *Dir;
    int ndir, tmp; // D = dimension, N = number of vectors
    OPT options;
    int i,j;
    double x;
    
    unsigned int seed;
    seed = (unsigned int) time(NULL);
    srand(seed);
    
    printf("ToDo\n - use natural sort\n - use stochastic gradient descent \n - add Gram-Schmidt normalization\n");

/* - - - - - - - -
 *  Get input data
 * - - - - - - - -
 */ 
    
    if (nrhs < 2)
    {
        mexErrMsgTxt("Sliced_Wasserstein_Projection: not enough inputs !\n->two arrays must be at least provided.");
    }
    X0Data = (mxArray *) prhs[0]; // input point cloud
    //X0 = mxGetPr(X0Data); // matrix
    options.D = mxGetM(X0Data); // number of rows
    options.N = mxGetN(X0Data); // number of columns
    
    YData  = (mxArray *) prhs[1]; // target point cloud
    if ( options.D != mxGetM(YData) || options.N != mxGetN(YData) ) {
        mexErrMsgTxt("Sliced_Wasserstein_Projection: incompatible arrays!\n->Input arrays must be compatible.");
    }
    //Y = mxGetPr(YData); // matrix
    
/* - - - - - - - -
 *  set options
 * - - - - - - - -
 */ 
    
    // a) verbose ?
    if (nrhs >= 8)
    {
        options.verbose = ( (int) * mxGetPr(prhs[7]) ) != 0; // 0 or 1
    }
    else 
    {
        options.verbose = 1; // default value
    }
    
    // b) set step descent parameter
    if (nrhs >= 3)
    {
        options.step = (double) * ( mxGetPr(prhs[2]) ); // gradient descent step (max = 1)
        if (options.step < 0) {
            options.step = 1;
            if (options.verbose)
                mexPrintf("Warning : step parameter should be positive ! set to : %d\n", options.step);
        }
    }
    else
    {
        options.step = 1; // default
        if (options.verbose)
            mexPrintf("Warning : Default step parameter : %g\n", options.step);
    }
    
    
    
    // c) Number of directions and direction set
    if (nrhs >= 4)
    {
        Dir  = mxGetPr(prhs[3]);
        ndir = mxGetN (prhs[3]);
        tmp  = mxGetM (prhs[3]);
        
        if ( tmp==1 && ndir==1 ) { // Dir is a scalar
            ndir = *Dir;
            if (options.verbose)
                mexPrintf("%d random directions generated ...\n", ndir);  
            
            // generer une matrice aléatoire et normalisée avec 'randn'
            Dir = (double *) mxCalloc(options.D * ndir, sizeof(*Dir));
            if (Dir==NULL)
                mexErrMsgTxt("Sliced_Wasserstein_Projection:memory allocation problem !\n");

            for (i=0;i<ndir;i++) { // columns
                for (j=0;j<options.D;j++) { // lines
                    Dir[j + i*options.D] = rand_normal(0,1); // Normal pdf
                    // Dir[j + i*options.D] = 2.0 * rand()/RAND_MAX -1.0; // Uniform pdf in [-1,1]^d
                }
            }
                    
            tmp = options.D;
        }
        
        if ( tmp != options.D )
            mexErrMsgTxt("Sliced_Wasserstein_Projection:array imcompatible !\n->Input array 'X0' and 'Dir' must have the same number of rows.");
            
        // normalization of direction set
        for (i=0;i<ndir;i++) { // column-wise loop
            x = 0.0;
            for (j=0;j<options.D;j++) { // rows = dimensions
                x += (Dir[j + i*(options.D)]) * (Dir[j + i*(options.D)]); // Dir[j,i]^2
            }
            // column normalization
            if (x>0) {
                for (j=0;j<options.D;j++) { // dimensions
                    Dir[j + i*(options.D)] /= sqrt(x); // L2 normalization
                }
            }
            else { // if null vector ...
                 Dir[i*(options.D)] = 1.0;
            }
        }
    }
    else
    {
        ndir = options.D; // default (minimal number of directions)
        if (options.verbose)
            mexPrintf("Warning : Default number of directions : %d ('Dir' is the identity matrix)\n", ndir);
        
        Dir = (double *) mxCalloc(options.D * ndir, sizeof(*Dir));
        if (Dir==NULL)
                mexErrMsgTxt("Sliced_Wasserstein_Projection:memory allocation problem !\n");

        for (i=0;i<ndir;i++) { // column-wise
            if (i<options.D)
                Dir[i + i*options.D] = 1.f;
            else
                Dir[i%(options.D) + i*options.D] = -1.f;
        }
    }
    
    // d) number of iterations
    if (nrhs >= 5)
    {
        options.niter = (int) * mxGetPr(prhs[4]);
        if (options.niter < 0) {
            options.niter = 100;
            if (options.verbose)
                mexPrintf("Warning : 'niter' parameter should be positive ! set to : %d\n", options.niter);
        }
    }
    else 
    {
        options.niter = 100; // default value
    }
    
    // e) other options
    if (nrhs >= 6)
    {
        options.flag_stochastic = (int) * mxGetPr(prhs[5]);
        if (options.flag_stochastic != 0) {
            options.flag_stochastic = 0;
            if (options.verbose)
                mexPrintf("Warning : 'flag_stochastic' parameter is not supported yet ! set to : %d\n", options.flag_stochastic);
        }
    }
    else 
    {
        options.flag_stochastic = 0; // default value
    }
    if (nrhs >= 7)
    {
        options.flag_hessian = (int) * mxGetPr(prhs[6]);
        if (options.flag_hessian != 0) {
            options.flag_hessian = 0;
            if (options.verbose)
                mexPrintf("Warning : 'flag_hessian' parameter is not supported yet ! set to : %d\n", options.flag_hessian);
        }
    }
    else 
    {
        options.flag_hessian = 0; // default value
    }
    
    // set option structure for algorithm
    options.ndir = ndir;
    options.Dir = Dir;
    options.flag_stochastic = 0;
    options.flag_hessian = 0;
    
    if (options.verbose) {
        mexPrintf("Options : \n");
        mexPrintf(" Size = %d x %d\n",options.D,options.N);
        mexPrintf(" niter = %d\n",options.niter);
        mexPrintf(" ndir = %d\n",options.ndir);
        mexPrintf(" step = %g\n",options.step);
        mexPrintf(" stochastic = %d\n",options.flag_stochastic);
        mexPrintf(" hessian = %d\n",options.flag_hessian);
        
        if (DEBUG) {
            mexPrintf(" Directions = \n");
            for (i=0;i<options.D;i++) { // lines
                for (j=0;j<options.ndir;j++) { // columns
                    mexPrintf("%g ",Dir[i+j*options.D]); // Dir[i,j]
                }
                mexPrintf("\n");
            }
        }
    }
    
    
    
/* - - - - - - - -
 *  Algorithm
 * - - - - - - - -
 */ 
    
    // initialize output(s) 'X'
    plhs[0] = mxDuplicateArray(prhs[0]);
    
    if (nlhs == 2)// optional output : energy
        plhs[1] = Sliced_Wasserstein_Projection(plhs[0],prhs[1],options); // (X,Y,options)
    else
        Sliced_Wasserstein_Projection(plhs[0],prhs[1],options);
        
    
    // mxFree(Dir); ONLY when using mexMakeMemoryPersistent 
    return;
}


// compute Sliced L2-Wasserstein Optimal transport of X onto Y; output is the Energy
mxArray * Sliced_Wasserstein_Projection(mxArray *X,const mxArray *Y,OPT options) {
    
    double *x, *y, * Dir = options.Dir;
    double step = options.step;
    int D = options.D, // D = dimension, 
        N = options.N, // N = number of vectors
        ndir = options.ndir, // Number of directions
        niter = options.niter, // Number of iterations
        flag_hessian = options.flag_hessian, // use hessian
        flag_stochastic = options.flag_stochastic; // use stochastic    
    int it,i,j,k,l;
    double s;
    
    mxArray * In[3];

    mxArray * Yproj = mxCreateDoubleMatrix(ndir, N, mxREAL);   
    double  * yproj = mxGetPr(Yproj);  
    
    mxArray * Ysort[2];
    double  * ysort; //= mxGetPr(Ysort[0]);
    
    mxArray * Xproj = mxCreateDoubleMatrix(ndir, N, mxREAL),
            * Xproj_presort = mxCreateDoubleMatrix(ndir, N, mxREAL); 
    double  * xproj = mxGetPr(Xproj),
            * xproj_presort = mxGetPr(Xproj_presort), // temporary array for sorting
            * x_tmp = mxCalloc(ndir * N, sizeof(* x_tmp)); // temporary array for permutation

    mxArray * Xsort[2] = {NULL,NULL};
    double  * xsort, // sorted value
            * xidx, // index for sorting
            * xidx_prev = mxCalloc(ndir * N, sizeof(* xidx_prev));
    int     * xidx_inv  = mxCalloc(ndir * N, sizeof(* xidx_inv ));
    
    mxArray * Grad = mxCreateDoubleMatrix(D, N, mxREAL); // gradient on sorted values
    double  * grad = mxGetPr(Grad);
    
    mxArray * E = mxCreateDoubleMatrix(niter, 1, mxREAL); 
    double  * e = mxGetPr(E);
    
    
    
    x = mxGetPr(X);
    y = mxGetPr(Y);
    
    // projection of the target N-point-cloud 'Y' on the ndir 'Dir' orientations set : Yproj = Dir' * Y
    for(j=0;j<N;j++) { // columns
        for(i=0;i<ndir;i++) { // rows
            yproj[i+j*ndir] = 0;
            for (k=0; k<D; k++) {
                yproj[i+(j*ndir)] += Dir[k + (i*D)] * y[k + (j*D)]; // yprok[i,j] = sum_k Dir[k,i] * Y[k,j]
            }
        }
    }
    
    In[0] = (mxArray *) Yproj; // array of values to be sorted
    In[1] = mxCreateDoubleScalar(2.0); // dimension along which is performed the sorting
    In[2] = mxCreateString("ascend"); // 'ascend' or 'descend' sorting 
    mexCallMATLAB(2, Ysort, 3, In, "sort"); // buit-in matlab sort is multi-threaded
    ysort = mxGetPr(Ysort[0]);
    
    // gradient descent
    if (options.verbose) mexPrintf("Gradient descent ...\n");
    for (it=0; it<options.niter; it++) {
        // count-down
        if (options.verbose && (20*it)%niter==0) { // 
            mexPrintf("[%2d%%]\n",(int) (100.0*it/(double)options.niter));
            mexEvalString("drawnow;"); // required to dump string
        }
        // step 1 : compute projection (scalar product between X[:,j] and Dir[:,i])
        if (DEBUG) mexPrintf("Step 1 ...\n");
        for(j=0;j<N;j++) { // columns
            for(i=0;i<ndir;i++) { // rows
                xproj[i+j*ndir] = 0.0;
                for (k=0; k<D; k++) {
                    xproj[i+(j*ndir)] += Dir[k + (i*D)] * x[k + (j*D)]; // yproj[i,j] = sum_k Dir[k,i] * Y[k,j]
                }
            }
        }
        
        // step 2 : sorting (using eventually acceleration using previous ordering)        
        if (DEBUG) mexPrintf("Step 2 ...\n");
        
        if (USE_PREVIOUS_ORDERING && it) { // sorting 'Xproj' according previous optimal permutation
            perform_permute_array(xproj, xproj_presort, xidx_prev, ndir, N);
            In[0] = (mxArray *) Xproj_presort; // array of values pre-sorted
        }
        else {
            In[0] = (mxArray *) Xproj; // array of values to be sorted
        }
        mxDestroyArray(Xsort[0]); // to prevent memory leak
        mxDestroyArray(Xsort[1]);
        mexCallMATLAB(2, Xsort, 3, In, "sort"); // buit-in matlab sort is multi-threaded
        xsort = mxGetPr(Xsort[0]);
        xidx  = mxGetPr(Xsort[1]); // NOTE: matlab index starts at 1 !
        // TEST : perform_permute_array(xproj, xsort, xidx, ndir, N);

        
        if (USE_PREVIOUS_ORDERING && it) {
            self_permute_array(xidx_prev, x_tmp, xidx, ndir, N);
            copy_array(xidx_prev, xidx, N * ndir);
            // the true sorting index is : xidx = xidx_prev o xidx
        }
        
        // step 2bis : compute inverted list 'xidx_inv', so that xidx_inv[i] = rank of xproj[i]
        // (corresponds to the matching Xproj[k,i]->Ysort[k,j] where j=xidx_inv[k+i*ndir])
        if (DEBUG) mexPrintf("Step 2.5 ...\n");
        for(j=0;j<N;j++) { // columns
            for(k=0;k<ndir;k++) { // rows
                xidx_inv[k + ((int) xidx[k + j*ndir]-1)*ndir] = j; 
                // corresponds to (f^{-1} o f = Id)
                // NOTE: '-1' is due to matlab index starting at '1' instead of '0'
            }
        }
        
        // step 3 : compute Sliced L2-Wasserstein gradient on sorted values
        if (DEBUG) mexPrintf("Step 3 ...\n");
        for(j=0;j<N;j++) { // columns : vectors
            for(i=0;i<D;i++) { // rows : dimension
                grad[i + j*D] = 0.0;
                for(k=0;k<ndir;k++) { // directions
                    l = xidx_inv[k+j*ndir]; // indices are now in [0,N-1]
                    grad[i + j*D] += (xproj[k+j*ndir] - ysort[k+l*ndir]) * Dir[i+k*D];
                }
                grad[i + j*D] /= (double) ndir; 
            }
        }
        
        // step 4 : perform gradient descent to update X
        if (DEBUG) mexPrintf("Step 4 ...\n");
        for(j=0;j<N;j++) { // columns : vectors
            for(i=0;i<D;i++) { // rows : dimension
                x[i+j*D] -= step * grad[i+j*D]; 
            }
        }
        
        // step 5 : compute Sliced L2-Wasserstein energy using 'xsort' and 'ysort'
        if (DEBUG) mexPrintf("Step 5 ...\n");
        e[it] = 0.0;
        for(j=0;j<N;j++) { // columns : vectors
            for(i=0;i<ndir;i++) { // rows : directions
                s = xsort[i+j*ndir] - ysort[i+j*ndir];
                e[it] += s*s;
            }
        }
        e[it] /= 2.0*ndir; 
        
        // update
        if (USE_PREVIOUS_ORDERING) {
            copy_array(xidx, xidx_prev, ndir * N);
        }
    }
    
    return E;
}

void copy_array(double *in,double *out,int i_max) {
    
    int i;
    for (i=0; i<i_max; i++) {
        out[i] = in[i];
    }
    
    return;
}

// permute from index array (cast as 'double') : Z = Z o p
void self_permute_array(double * Z, double * z, double * p, int ndir, int N) {
    int i,j,k;
    copy_array (Z,z,ndir*N);
    
    for(j=0;j<N;j++) { // col
        for(i=0;i<ndir;i++) { // rows
            k = (int) p[i+j*ndir] - 1;
            Z[i + j*ndir] = z[i + k*ndir];
        }
    }
    
    return;
}
// permute : out = in o p
void perform_permute_array(double * in, double * out, double * p, int ndir, int N) {
    int i,j,k;
    
    for(j=0;j<N;j++) { // col
        for(i=0;i<ndir;i++) { // rows
            k = (int) p[i+j*ndir] - 1;
            out[i + j*ndir] = in[ i + k*ndir ];
        }
    }
    
    return;
}

// Box-Muller function to generate normal random value
double rand_normal(double mean, double stddev) {
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached) {
        // choose a point x,y in the unit circle uniformly at random
        double x, y, r;
        do {
            // scale two random integers to doubles between -1 and 1
            x = 2.0 * rand() / RAND_MAX - 1;
            y = 2.0 * rand() / RAND_MAX - 1;
            r = x*x + y*y;
        } while (r == 0.0 || r > 1.0);
        
        
        // Apply Box-Muller transform on x, y
        double d = sqrt(-2.0*log(r)/r);
        double n1 = x*d;
        n2 = y*d;
        // scale and translate to get desired mean and standard deviation
        double result = n1*stddev + mean;
        n2_cached = 1;
        return result;
        
    }
    else {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}
 
