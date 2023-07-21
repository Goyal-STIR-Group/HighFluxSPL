%% Simulating electronics and detector dead time
% Joshua Rapp
% Boston University
%
% Comparison of detection PDF accounting for only detector dead time (Td)
% or both detector and electronics dead times (Td and Te)
%
% For code with Td and Te, please cite
% J. Rapp, Y. Ma, R. Dawson, and V. Goyal, 
% “High-flux single-photon lidar,” Optica, vol. 8, no. 1, 
% pp. 30–39, Nov. 2021, doi: 10.1364/optica.403190.
%
% For code with Td only, please cite
% J. Rapp, Y. Ma, R. M. A. Dawson, and V. K. Goyal, 
% “Dead time compensation for high-flux ranging,” 
% IEEE Transactions on Signal Processing, pp. 1–1, Oct. 2019, 
% doi: 10.1109/TSP.2019.2914891.


clear; close all; clc;
addpath('utils');

%% Generate Fig. 4 in “High-flux single-photon lidar”
% Parameters
Nr = 5e5;       % Number of illumination periods
S = 1;          % Mean # signal photons per period
B = 5;          % Mean # background photons per period

% Make sure Td < Te <= Tr <= 2Td
Tr = 100;       % Illumination rep. period (ns)
Td = 55;        % Detector dead time (ns)
Te = 80;        % Electronics dead time (ns)
Tbin = 0.05;    % Bin size [ns]

Ztrue = 75;     % True signal position (ns)
sigP = 0.5;     % Pulse width

% Compute attenuated detection sequence from same photon arrivals?
LFflag = true;

% Unused parameters - to add uncertainty to detection
Te_subMean = 0;     % mean of exponential distribution (ns)
Td_sig = 0;         % std of Gaussian distribution (ns)
detDecay = 0;       % mean of exponential 

% Perform photon detection sequence simulation
[ LFDetTimes,HFDetTimes ] = ...
    fcn_gen_DE_Dets3( Nr,S,B,Tr,Td,Ztrue,sigP,Te,Te_subMean,Td_sig,detDecay,LFflag);

% Histogram discretization
bins = Tbin/2 : Tbin : Tr-Tbin/2; % bin centers [ns]

% True detection time PDF
lam_true =S*normpdf(bins,Ztrue,sigP)+B/Tr;
f_arriv_true = lam_true/sum(lam_true*Tbin);

% PDF with electronics and detector dead times
f_TdTe =  fcn_PDF_PileupTdTe(Tbin, lam_true, Tr, Td, Te);
f_det_TdTe = f_TdTe/sum(f_TdTe*Tbin);

% PDF with detector dead time only
f_Td =  fcn_PDF_PileupTd(lam_true, Td, Tr, Tbin);
f_det_Td = f_Td/sum(f_Td*Tbin);

% Plot
cmap = lines(5);
figpos = [10 10 8 6.8];
hfigHF = figure('units','centimeters','position',figpos); 
histogram(mod(HFDetTimes,Tr),50,'Normalization','pdf'); 
hold on; set(gca,'yscale','log','fontsize',12);
plot(bins,f_arriv_true,'linewidth',2,'linestyle','-','color',cmap(2,:));
plot(bins,f_det_TdTe,'k-','linewidth',2);
plot(bins,f_det_Td,'linewidth',2,'linestyle','--','color',cmap(4,:));
xlabel('Time [ns]');
