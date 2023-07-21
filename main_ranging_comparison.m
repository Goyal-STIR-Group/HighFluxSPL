% HF PDF Estimation
% Joshua Rapp
% Boston University
% July 01, 2019

clear; clc; close all;
savePlots = false;

acqDate = '2019_01_28';

% Choose data set index (2-14)
dSIndx = 2;
% Choose downsampling rate (to speed up processing)
downSampFactor = 4;

% File management
addpath('utils');
pname = ['Data/' acqDate '/'];
HFacqNum = 2:14;
LFacqNum = ones(13,1);
laserIntensity = 100*ones(13,1);

T_ho = [48;48;57;66;75;78;81;87;96;99;120;147;180;198];
Td_ests = T_ho(2:end)+2;%[50;58.5;68;76;79;82.5;89;98;101;122;149;181;198.5];
Bhf_ests = [1.87;1.79;1.75;1.65;1.65;1.72;1.72;1.68;1.65;1.63;1.6;1.58;1.57];

fBaseHF = ['FGS_td' sprintf('%03i',T_ho(HFacqNum(dSIndx))) '_' acqDate '_acq'];
fBaseLF = ['FGS_td' sprintf('%03i',T_ho(LFacqNum(dSIndx))) '_' acqDate '_acq'];
disp(['Results for acquisition with ' num2str(T_ho(HFacqNum(dSIndx))) ' ns hold-off time']);

Te = 80; % electronics dead time [ns]

%% Load High Flux Experimental Data

% load data
load([pname fBaseHF sprintf('%i',HFacqNum(dSIndx)) '.mat']);

binResHF = binRes;
binResDS = binResHF*downSampFactor;
numBins = floor(tr/binResDS);
binEdges = ((0:numBins))*binResDS*1e9;
binCenters = ((0:numBins-1)+0.5)*binResDS*1e9;
Tbin = binResDS*1e9;

detTimeHF = detTime;%*binResHF*1e9;
frameNumHF = frameNum;
OD_HF = OD;
Tr = tr*1e9;

detHistHF = histcounts(detTimeHF,binEdges);
% figure; semilogy(binCenters,detHistHF); title(['T_{hold-off} = ' num2str(T_holdoff) ' ns']);
% ylabel('Detection Count'); xlabel('Time [ns]');

% Estimate total flux Lambda
absDetTimeHF = frameNumHF*Tr+detTimeHF;
Td = T_holdoff+2;%min(T_holdoff+2,min(diff(absDetTimeHF)));   % Detector dead time is approx. 2 ns longer than hold off

if Td<=Te
    Tmin = Te+Td;
else
    Tmin = Td;
end
diffAbsDetTimeHFall = diff(absDetTimeHF);
diffAbsDetTimeHF = diffAbsDetTimeHFall(diffAbsDetTimeHFall>Tmin);
Rhf = floor((diffAbsDetTimeHF-Tmin)/Tr);
LambdaEstHF = -log(sum(Rhf)/(length(Rhf) + sum(Rhf)));

interPeriodsHF = 0:max(Rhf);
HFdecay = (1-exp(-LambdaEstHF))*exp(-LambdaEstHF*interPeriodsHF);

hfEdges = 0:3:4.5*Tr;

% Acquisition statistics
totCountHF = length(detTimeHF);
disp(['HF counts = ' num2str(totCountHF)]);

figpos = [10 10 10 6.6];
hfigIT = figure('units','centimeters','position',figpos);
histogram(diffAbsDetTimeHFall,hfEdges); title('HF Interdetection times');
xlabel('Time [ns]'); ylabel('Counts'); set(gca,'yscale','log','fontsize',14);

hfigIP = figure('units','centimeters','position',figpos+[figpos(3) zeros(1,3)]);
histogram(Rhf,'normalization','pdf');
hold on; plot(interPeriodsHF,HFdecay,'linewidth',1.5);
title('HF Interdetection Periods'); xlabel('# Periods'); ylabel('Frequency');
set(gca,'yscale','log','fontsize',14);
legend('Measured','Estimated');


%% Get rid of large data (no longer needed)
clear Rhf absDetTimeHF detTimeHF diffAbsDetTimeHF diffAbsDetTimeHFall frameNumHF

%% Load Low Flux Experimental Data
%lam =S*normpdf(bins,Ztrue,sigP)+B/Tr;
load([pname fBaseLF sprintf('%i',LFacqNum(dSIndx)) '.mat']);

binResLF = binRes;
detTimeLF = detTime;%*binResLF*1e9;
frameNumLF = frameNum;
OD_LF = OD;
OD_Diff = OD_LF-OD_HF;

detHistLF = histcounts(detTimeLF,binEdges);

% Lambda Estimation from Low-Flux measurement
absDetTimeLF = frameNumLF*Tr+detTimeLF;
diffAbsDetTimeLFall = diff(absDetTimeLF);
diffAbsDetTimeLF = diffAbsDetTimeLFall(diffAbsDetTimeLFall>Tmin);
Rlf = floor((diffAbsDetTimeLF-Tmin)/Tr);
LambdaEstLF = (-log(sum(Rlf)/(length(Rlf) + sum(Rlf))));%*10^OD_Diff;

interPeriodsLF = 1:4:1000;
LFdecay = (1-exp(-LambdaEstLF))*exp(-LambdaEstLF*interPeriodsLF);

lfEdges = 0:3:4.5*Tr;

%% lambda(t) estimation
% Load pulse shape calibration + normalize
load([pname 'laser_calib_' sprintf('%i',laserIntensity(dSIndx)) '.mat']);
binResPS = binRes;

downSampFactorPS = downSampFactor*binResHF/binResPS;
pulseShape = (smoothPulseShape-min(smoothPulseShape))/sum(smoothPulseShape-min(smoothPulseShape));
numPulseDS = floor(numel(pulseShape)/downSampFactor);
pulseShapeDS = zeros(numPulseDS,1);
for ii = 1:downSampFactor
    dsTemp = downsample(pulseShape,downSampFactor,ii-1);
    pulseShapeDS = pulseShapeDS + dsTemp(1:numPulseDS);
end
pulseWindow = 90:90+numPulseDS-1;

TdEst = Td_ests(dSIndx);
TeEst = Te;
empPDF = detHistHF/sum(detHistHF);
empCDF = cumsum(empPDF);

numTest = 11;
BhfTests = Bhf_ests(dSIndx)+linspace(-0.05,0.05,numTest);
KSds = zeros(numTest,1);
piTHats = zeros(numBins,numTest);

disp('Optimizing Kolmogorov-Smirnov statistic:');
parfor jj = 1:numTest
    disp(['... ' num2str(jj) ' of ' num2str(numTest)]);
    BhfTEst = BhfTests(jj);%Bhf_ests(dSIndx);%BEstHF;%
    ShfTEst = LambdaEstHF-BhfTEst;
    LhfEst = BhfTEst+ShfTEst;
    
    intArriv1 = BhfTEst*ones(1,numBins)/Tr;
    intArriv1(pulseWindow) = ShfTEst*pulseShapeDS/Tbin+BhfTEst/Tr;
%     lambdaLF = detHistLF/sum(detHistLF*Tbin)*LhfEst;
    
    % Compute transition probability matrix
    lam = intArriv1;
    
    if TdEst < TeEst
        piTHat =  fcn_PDF_PileupTdTe(Tbin, lam, Tr, TdEst, TeEst);
    else
        piTHat = fcn_PDF_PileupTd(lam,TdEst,Tr,Tbin);
    end
    piTHats(:,jj) = piTHat;
    
    % Compute Closeness of fit statistics
    predCDF = cumsum(piTHat)*Tbin;
    
    % Kolmogorov Smirnov statistic
    KSds(jj) = max(abs(empCDF(:)-predCDF(:)));
end
disp('... done.');

% KS Optimization results
[~, minDex] = min(KSds);
BhfEst = BhfTests(minDex);
ShfEst = LambdaEstHF-BhfEst;
piHat = piTHats(:,minDex);

disp(['Estimated B = ' num2str(BhfEst)]);
disp(['Estimated S = ' num2str(ShfEst)]);

% Predicted PDFs
intArrivEst = BhfEst*ones(1,numBins)/Tr;
intArrivEst(pulseWindow) = ShfEst*pulseShapeDS/Tbin+BhfEst/Tr;

piHatTdOnly = fcn_PDF_PileupTd(intArrivEst,TdEst,Tr,Tbin);%
piHatTdOnly = piHatTdOnly/sum(piHatTdOnly*Tbin);

piHatTeOnly = fcn_PDF_PileupTd(intArrivEst,TeEst,Tr,Tbin);%
piHatTeOnly = piHatTeOnly/sum(piHatTeOnly*Tbin);

predCDF = cumsum(piHat)*Tbin;


%% Plots 
% Low Flux
cmap = lines(5);
insetXMin = 1.0e-9;
insetXMax = 6.0e-9;
insetBins = round(insetXMin/binResDS):round(insetXMax/binResDS);

figpos = [10 10 10 7.5];
hfigHF = figure('units','centimeters','position',figpos);
semilogy(binCenters,detHistLF/sum(detHistLF*Tbin),'linewidth',1.5);
hold on; semilogy(binCenters,intArrivEst/sum(intArrivEst*Tbin),'linewidth',1.5); 
plot(insetXMin*1e9*ones(2,1),[1e-3;1],'--','color',cmap(4,:),'linewidth',1.1);
plot(insetXMax*1e9*ones(2,1),[1e-3;1],'--','color',cmap(4,:),'linewidth',1.1);
xlabel('Time [ns]'); 
ylim([5*10^-3 1]); 
set(gca,'fontsize',12);
legend('Measured histogram','Smoothed estimate','location','northeast')

axes('position', [0.4 0.38 0.4 0.32]);
semilogy(binCenters(insetBins),detHistLF(insetBins)/sum(detHistLF*Tbin),'linewidth',1.5);
hold on; semilogy(binCenters(insetBins),intArrivEst(insetBins)/sum(intArrivEst*Tbin),'linewidth',1.5); 
xlim([insetXMin insetXMax]*1e9); ylim([6e-3 1]); 
set(gca,'fontsize',10);

% High Flux
insetYMin = 6e-3;
insetYMax = 0.013;

figpos = [10 10 10 7.5];
hfigHF = figure('units','centimeters','position',figpos);
semilogy(binCenters,detHistHF/sum(detHistHF*Tbin),'linewidth',1.5);
hold on; 
% semilogy(binCenters,piHatTdOnly,':','color',cmap(5,:),'linewidth',1.8); 
% semilogy(binCenters,piHatTeOnly,'k-.','linewidth',1.5);
semilogy(binCenters,piHat,'color',cmap(2,:),'linewidth',1.5);
xlabel('Time [ns]'); 
ylim([5*10^-3 1]); %axis tight;
set(gca,'fontsize',12);

axes('position', [0.3 0.43 0.5 0.47]);
semilogy(binCenters,detHistHF/sum(detHistHF*Tbin),'linewidth',1.5);
hold on; 
% semilogy(binCenters,piHatTdOnly,':','color',cmap(5,:),'linewidth',1.8); 
% semilogy(binCenters,piHatTeOnly,'k-.','linewidth',1.5);
semilogy(binCenters,piHat,'color',cmap(2,:),'linewidth',2);

ylim([insetYMin insetYMax]); yticks([1e-2]); 
ax = gca;
ax.YAxis.Exponent = 0;
set(gca,'fontsize',10);
