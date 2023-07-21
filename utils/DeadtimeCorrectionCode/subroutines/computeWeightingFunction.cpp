// % Authors: Simon Christoph Stein and Sebastian Isbaner
// % Year: 2016
// % E-Mail: scstein@phys.uni-goettingen.de

// mex
#include "mex.h"

// stdlib
#include <math.h>
#include <vector>
#include <iostream>

// our headers
#include "mexUtil.h"

using namespace std;


// -- Type definitions --
typedef Array1D_t<double> Array1D;

// -- Prototypes -- //

// --  Global Variables  -- //
mxArray** plhs;  // Access to left hand side is global, so we can output to it in functions


// Function definitions 

// Access array 'a' assuming periodicity of 'a' for indices that are negative or larger then the size of 'a'.
// E.g. pAcc(a,-1) = a[a.nElements-1];
double pAcc(Array1D& a, int idx);

// Integrate array in the interval [lowerBound,upperBound) assuming a is periodic (see pAcc above);
double intArr(Array1D& a, int lowerBound, int upperBound);

// Simplest possible algorithm to compute w(t)
void computeW_simple(Array1D& w, Array1D& h,Array1D& k, int E,int D, double epsilonP);

// Fast computation
// Compute w(t) using running sums to evaluate outer integrals
void computeW_fast(Array1D& w, Array1D& h,Array1D& k, int E,int D, double epsilonP);

// Fastest computation
// Compute w(t) using running sums to evaluate outer and inner integrals
// and use the periodicity of k to replace integrals over one period by epsilonP
void computeW_fastest(Array1D& w, Array1D& h,Array1D& k, int E,int D, double epsilonP);

/*
% Usage: [ w ] = computeWeightingFunction(h,k,E,D,epsilon)
%  Computes the weighting function w(t) which can be used to recover the
%  true photon hit rate k(t) from the measured curve h(t) = w(t)*k(t)
%  taking dead-time effects of electronics and detector into account.
%
%  To do this, intialize k(t) = h(t),
%  -->
%  | (x) compute w(t)
%  | (x) k(t)=h(t)/w(t)
%  -- iterate
%
% For more info, see 
% "Dead-time correction of fluorescence lifetime measurements", 
% Sebastian Isbaner, Narain Karedla et al. (submitted)
%
% Input:
%  h(t) - 1D double array, the measured decay curve.
%  k(t) - 1D double array, current estimate of true photon hit rate. 
%         Must be normalized to Int_0^P k(t) dt = epsilonP, with P period of k(t)!
%  E - Integer. Electronics dead-time.
%  D - Integer. Detector dead-time.
%  epsilonP - Double. Average number of hitting photons per excitation  cycle.
%
% Authors: Simon Christoph Stein and Sebastian Isbaner
% Year: 2016
% E-Mail: scstein@phys.uni-goettingen.de
*/
void mexFunction(int nlhs, mxArray* plhs_arg[], int nrhs, const mxArray* prhs[])
{    
    // Access to left hand side is global, so we can output to it in functions
    plhs = plhs_arg;
   
    /* Check for proper number of arguments. */
    if(nrhs<5) {
        mexErrMsgTxt("Five inputs required w(h,k,E,D,epsilon).");
    }    
    
    /* Map access to input data */
    Array1D h( prhs[0] );
    Array1D k( prhs[1] );
    int E = int(*mxGetPr(prhs[2]) + 0.5);
    int D = int(*mxGetPr(prhs[3]) + 0.5);
    double epsilonP = *mxGetPr(prhs[4]);
    
    int mode = 1; // default mode
    if(nrhs==6)
    {    
        mode = int(*mxGetPr(prhs[5]) + 0.5);
    }
    
    // -- Reserve Matlab output memory -- //
    mwSize dim_out[1] = { h.nElements };
    plhs[0] = mxCreateNumericArray( 1, dim_out , mxDOUBLE_CLASS, mxREAL);
    
    Array1D w ( plhs[0] );
    
    // Compute w(t)
    switch(mode)
    {
        case 1:
            computeW_fastest(w, h,k,E,D,epsilonP); break;
        case 2:
            computeW_fast(w, h,k,E,D,epsilonP); break;
        case 3:
            computeW_simple(w, h,k,E,D,epsilonP); break;
        default:
            mexWarnMsgTxt("Unknown mode, defaulting to mode=1 (fastest algorithm).");
            computeW_fastest(w, h,k,E,D,epsilonP);
            break;
    }    
}


// Access array 'a' assuming periodicity of 'a' for indices that are negative or larger then the size of 'a'.
// E.g. pAcc(a,-1) = a[a.nElements-1];
double pAcc(Array1D& a, int idx)
{
  const int Period = a.nElements;
  int remainder = idx%Period;
  
  if (remainder==0)
      return a[0];
  if (remainder>0)
      return a[remainder];
  if (remainder<0)
      return a[Period  + remainder];
}


// Integrate array in the interval [lowerBound,upperBound) assuming a is periodic (see pAcc above);
double intArr(Array1D& a, int lowerBound, int upperBound)
{
    double sum = 0;
    if (lowerBound < upperBound)
        for(int idx = lowerBound; idx<upperBound; ++idx)
        {
            sum += pAcc(a,idx);
        }
    else
    {
        for(int idx = upperBound; idx<lowerBound; ++idx)
        {
            sum -= pAcc(a,idx);
        }
    }
    return sum;
}

// Simplest possible algorithm to compute w(t)
void computeW_simple(Array1D& w, Array1D& h,Array1D& k, int E,int D, double epsilonP)
{
    int P = h.nElements;
    
    // For every point in time
    for (int t = 0; t<P; ++t)
    {
        double first_sum = exp(-intArr(k,t-D,t)) * intArr(h,t-E-D,t-E);

        double second_sum = 0;
        for (int tp = t-E-D-P; tp<t-E-D; ++tp)
        {
            second_sum += pAcc(h,tp) * exp(-intArr(k,tp+E,t));
        }

        w[t] = first_sum + 1./(1.-exp(-epsilonP)) * second_sum;
    }
}

// Fast computation
// Compute w(t) using running sums to evaluate outer integrals
void computeW_fast(Array1D& w, Array1D& h,Array1D& k, int E,int D, double epsilonP)
{
    int P = h.nElements;
    
    // Initialize running sums
    double int_tD_t_k = intArr(k,0-D,0);   
    double int_tED_tE_h = intArr(h,0-E-D,0-E);
    
    double first_sum = exp(-int_tD_t_k) * int_tED_tE_h;
    
    
    double second_sum = 0;
    for (int tp = 0-E-D-P; tp<0-E-D; ++tp)
    {
        second_sum += pAcc(h,tp) * exp(-intArr(k,tp+E,0));
    }
    
    w[0] = first_sum + 1./(1.-exp(-epsilonP)) * second_sum;
    
    // For every point in time
    for (int t = 1; t<P; ++t)
    {
        // Update running sums
        int_tD_t_k = int_tD_t_k - pAcc(k,(t-1)-D) + pAcc(k,(t-1));
        int_tED_tE_h = int_tED_tE_h - pAcc(h,(t-1)-E-D) + pAcc(h,(t-1)-E);
        
        first_sum = exp(-int_tD_t_k) * int_tED_tE_h;
        
        second_sum -= pAcc(h,(t-1)-E-D-P) * exp(-intArr(k,(t-1)-D-P,(t-1)));
        second_sum += pAcc(h,(t-1)-E-D) * exp(-intArr(k,(t-1)-D,t));
        
        // Int_{t-E-D-P}^{t-E-D} dt' h(t') exp(-Int_{t'+E}^t dt'' k(t''))
        // t->(t+1) --> Int_{(t+1)-E-D-P}^{(t+1)-E-D} dt' h(t') exp(-Int_{t'+E}^(t+1) dt'' k(t''))
        //            = [Int_{(t+1)-E-D-P}^{(t+1)-E-D} dt' h(t') exp(-Int_{t'+E}^t) dt'' k(t''))] * exp(-Int_t^{t+1} dt'' k(t''))
        // The [..] part can be updated like other running sums and must be multiplied by exp(-Int_t^{t+1} dt'' k(t''))
        second_sum *= exp(-pAcc(k,(t-1))); 
                
        w[t] = first_sum + 1./(1.-exp(-epsilonP)) * second_sum;
    }
}


// Fastest computation
// Compute w(t) using running sums to evaluate outer and inner integrals
// and use the periodicity of k to replace integrals over one period by epsilonP
void computeW_fastest(Array1D& w, Array1D& h,Array1D& k, int E,int D, double epsilonP)
{
     int P = h.nElements;
    
    // Initialize running sums
    double int_tD_t_k = epsilonP * D/P + intArr(k,0-(D%P),0); // Check if D>P, use periodicity of k and Int_0^P dt k(t) = epsilonP
    double int_tED_tE_h = intArr(h,0-E-D,0-E);
    
    double first_sum = exp(-int_tD_t_k) * int_tED_tE_h;
    
    
    double second_sum = 0;
    double int_tpE_0_k = epsilonP*(1+D/P) + intArr(k,0-(D%P),0); // Note: Equal to intArr(k,0-D-P,0) using periodicity of k and Int_0^P dt k(t) = epsilonP
    for (int tp = 0-E-D-P; tp<0-E-D; ++tp)
    {
        second_sum += pAcc(h,tp) * exp(-int_tpE_0_k);
        int_tpE_0_k -= pAcc(k,tp+E);
    }
    
    w[0] = first_sum + 1./(1.-exp(-epsilonP)) * second_sum;
    
    // Initialize intArr(k,(t-1)-D,t) with t=1
    // Check if D>P, use periodicity of k and Int_0^P dt k(t) = epsilonP
//     double int_t1D_t_k = epsilonP * floor(D/P) + intArr(k,(1-1)-(D%P),1); 
    double int_t1D_t_k = int_tD_t_k + pAcc(k,0); 
    
    // For every point in time
    for (int t = 1; t<P; ++t)
    {
        // Update running sums
        int_tD_t_k = int_tD_t_k - pAcc(k,(t-1)-D) + pAcc(k,(t-1));
        int_tED_tE_h = int_tED_tE_h - pAcc(h,(t-1)-E-D) + pAcc(h,(t-1)-E);
        
        first_sum = exp(-int_tD_t_k) * int_tED_tE_h;
        
      
        // Note: We split Int_{t-1-D-P}^{t-1} dt' k(t') 
        //                = Int_{t-1-D-P}^{t-1-D} dt' k(t') + Int_{t-1-D}^{t-1} dt' k(t')
        //                =           epsilonP              + Int_{t-1-D}^{t} dt' k(t') - Int_{t-1}^t dt' k(t')
        // This way we only need to update one running sum for  Int_{t-1-D}^{t} dt' k(t').
        // BEFORE
//         second_sum -= pAcc(h,(t-1)-E-D-P) * exp( -intArr(k,(t-1)-D-P,(t-1))  );
//         second_sum += pAcc(h,(t-1)-E-D) * exp(-intArr(k,(t-1)-D,t)); 
        // AFTER
        second_sum -= pAcc(h,(t-1)-E-D-P) * exp( - (int_t1D_t_k + epsilonP - pAcc(k,t-1))  );
        second_sum += pAcc(h,(t-1)-E-D) * exp(-int_t1D_t_k);
        
        // Int_{t-E-D-P}^{t-E-D} dt' h(t') exp(-Int_{t'+E}^t dt'' k(t''))
        // t->(t+1) --> Int_{(t+1)-E-D-P}^{(t+1)-E-D} dt' h(t') exp(-Int_{t'+E}^(t+1) dt'' k(t''))
        //            = [Int_{(t+1)-E-D-P}^{(t+1)-E-D} dt' h(t') exp(-Int_{t'+E}^t) dt'' k(t''))] * exp(-Int_t^{t+1} dt'' k(t''))
        // The [..] part can be updated like other running sums and must be multiplied by exp(-Int_t^{t+1} dt'' k(t''))
        second_sum *= exp(-pAcc(k,(t-1)));
        
        int_t1D_t_k = int_t1D_t_k + pAcc(k,t) - pAcc(k,(t-1)-D);
                
        w[t] = first_sum + 1./(1.-exp(-epsilonP)) * second_sum;
    }
}



