function [tcspc_corrected,tcspc_standard,epsilonP] = TCSPC_DTcorrected(tPhoton,syncPeriod,DT_electronics,DT_detector,flag)
% Program for calculating both the "standard" and the unbiased, dead-time
% corrected TCSPC curve from measured photon arrival times. Crucial
% parameters are the dead-times of the detector and the electronics, as
% well as the sync period.
% All inputs have to be in the same time units, namely the width of one TCSPC-
% bin (typically tens of picoseconds). Since they will be used as indices to
% access matrix elements, they have to be whole numbers (although they may be
% of the data type "double").
% For example, with a sync period of 50ns, an electronics deadtime of 80ns,
% a detector deadtime of 33ns and TCSPC bins of 0.064ns, we have:
% syncPeriod=782; DT_electronics=1250; DT_detector=516; (always round up)
%
% Please note that the origin of our TCSPC-histograms is chosen such that
% photons with arrival times that correspond to multiples of the sync period
% will be sorted in the first bin of the resulting TCSPC-histograms.
%
% If flag is set to 1, the folder containing the subroutines has to be
% added to the MATLAB path manually or the subroutines have to be copied
% into the same folder as TCSPC_DTcorrected.m.
%
% Input:
% tPhoton        - absolute arrival times of the photons (vector)
% syncPeriod     - duration of one sync cycle            (scalar)
% DT_electronics - dead time of the electronics          (scalar)
% DT_detector    - dead time of the detector             (scalar)
% flag           - should the subfolder be added manually? (scalar)
%
% Output:
% tcspc_corrected - corrected TCSPC-curve; normalized such that
%                   sum(tcspc_corrected) = epsilon*P = "the average 
%                   number of photons per sync cycle"
% tcspc_standard  - TCSPC-curve obtained by binning without correction; not
%                   normalized, i.e. sum(tcspc_standard) = numel(tPhoton)
% epsilonP        - average number of photons emitted per excitation period
%                   (can be used to get a corrected intensity value for the
%                    pixel via intensity = epsilonP * number of sync cycles)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add subfolder 'subroutines' to the Matlab path (for running session only)
if nargin<5 || isempty(flag) || flag~=1
    fullPathToThisFile = mfilename('fullpath'); % Path to m-File executing this function
    [path,~,~] = fileparts(fullPathToThisFile);
    addpath(genpath([path,filesep,'subroutines'])); 
end

% get the "standard" TCSPC-curve by finding the TCSPC-bin of each photon  %
% and histogramming                                                       %
tcspc = mod(tPhoton,syncPeriod)+1; % we always start at time 1 instead of 0
tcspc_standard = mHist(tcspc,1:syncPeriod);

% estimate the true average number of photons per sync period, epsilon*P, %
% from the inter-photon time-differences                                  %
iptd = diff(tPhoton); % inter-photon time-difference
N = (1:10)'; mN = zeros(size(N));
rangeLeft = DT_electronics+DT_detector;
rangeRight= rangeLeft+syncPeriod;
for i = 1:numel(N) % integrate iptd over N intervals of width "syncPeriod"
    mN(i) = numel(find(iptd>rangeLeft & iptd<=rangeRight));
    rangeLeft = rangeLeft +syncPeriod;
    rangeRight= rangeRight+syncPeriod;
end
if any(mN==0)
    mN=mN(1:find(mN==0,1)-1);
    N = N(1:numel(mN));
end
weight = mN; % weighting vector, typically 1/error^2
NN=[N,ones(size(N))]; % model matrix with ones for constant offset
b = (NN'*diag(weight)*NN)\NN'*diag(weight)*(-log(mN));
epsilonP = b(1);

% average number of DETECTED & COUNTED photons per excitation period
epsilonP_apparent = numel(tPhoton)/(max(tPhoton)-min(tPhoton))*syncPeriod;
% the number of detected & counted photons can never be larger than the
% true number of emitted photons
if epsilonP<epsilonP_apparent
    epsilonP = epsilonP_apparent;
end

% get the dead-time corrected TCSPC curve
tcspc_corrected = DeadTimeCorrection(tcspc_standard,DT_electronics,DT_detector,epsilonP);
end
