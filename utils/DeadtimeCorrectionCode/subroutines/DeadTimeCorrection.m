function meas_corr = DeadTimeCorrection(meas,DT_electronics,DT_detector,epsilonP)
% Program for recovering an unbiased decay-curve from a measured TCSPC curve

% Input:
%    meas          - measured decay curve
%   DT_electronics - electroncis dead-time in time.units of TCSPC channels
%   DT_detector    - detector dead-time in time.units of TCSPC channels
%   epsilonP       - average number of hitting photons per excitation cycle
%
% Output:
%   meas_corr - the corrected photon hit rate normalized such that sum(back)=epsilonP
%

% If the maximum relative error between two iterations drops below 1-e3, we stop iterating.
REL_ERR_THRESH = 1e-3; 
% Maximum number of iterations to perform
MAX_IT = 20; 

meas = double(meas(:))/sum(meas)*epsilonP; % Take input as row vector and normalize
meas_corr = meas; % k(t). This must be normalized to sum(k(t))==epsilonP! (important!)
last_estimate = meas_corr;

for iteration=1:MAX_IT
    w = computeWeightingFunction(meas,meas_corr,DT_electronics,DT_detector,epsilonP);
    meas_corr = meas./w;
    meas_corr = meas_corr/sum(meas_corr)*epsilonP; % Normalize
    
    % Convergence check
    % If the maximum relative error drops below 1-e3, we stop iterating.
    relative_error = max(abs(meas_corr-last_estimate)./last_estimate);
    if(relative_error< REL_ERR_THRESH)
        return
    end
    last_estimate = meas_corr;
end
warning('MAX_ITERATION reached.')
end



