function [ w ] = computeWeightingFunction(h,k,E,D,epsilonP, mode)
% Usage: [ w ] = computeWeightingFunction(h,k,E,D,epsilon, mode)
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
%  h(t) - 1D double array, the measured decay curve
%  k(t) - 1D double array, current estimate of true photon hit rate. 
%         Must be normalized to Int_0^P k(t) dt = epsilonP, with P period of k(t)!
%  E - Integer. Electronics dead-time.
%  D - Integer. Detector dead-time.
%  epsilonP - Double. Average number of hitting photons per excitation  cycle.
%  mode (Optional) - Integer. Algorithm used for computing w(t). For debugging only, defaults to mode=1 if not given.
%         mode=1: Fastest algorithm
%         mode=2: Fast algorithm
%         mode=3: Slow algorithm
%      This is for testing if the faster (and more cryptic) algorithms do their job properly.
%
% Authors: Simon Christoph Stein and Sebastian Isbaner
% Year: 2016
% E-Mail: scstein@phys.uni-goettingen.de
end

