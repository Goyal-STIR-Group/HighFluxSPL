%% Experimental Dead Time Imaging
% Joshua Rapp
% Boston University
% July 10, 2019
%
% Markov Chain-based PDF prediction accounting for both detector and 
% electronics dead times (Td and Te).
%
% Generates Fig. 7 in the paper.
%
% Please cite
% J. Rapp, Y. Ma, R. Dawson, and V. Goyal, 
% “High-flux single-photon lidar,” Optica, vol. 8, no. 1, 
% pp. 30–39, Nov. 2021, doi: 10.1364/optica.403190.
%
% Processing acquisitions collected on May 24, 2019 at MIT.

clear; close all; clc;

%% Acquisition details
acqDate = '2019_05_24';
pname = 'Data/';
blAcq = 2;
lfAcq = 2;
bgndAcq = 4;
hfAcq = 6;

% Choose downsampling rate
downSampFactor = 2;

% Figure save settings
printAllFigs = 0;
dtime = datestr(datetime,30);
prname = ['../Figures/' acqDate '/'];
addpath('utils');

%% Load Baseline Data
BgndBL = 0.01;  % Guess of background in baseline data

load([pname acqDate '_dt_data_acq' num2str(blAcq) '_Histograms.mat']);
SigBL = LambdaMat;
[numY, numX] = size(LambdaMat);

load([pname acqDate '_dt_data_acq' num2str(blAcq) '_Details.mat'],'Resolution','SyncRate');

% Downsample time dimension for faster processing
% histMatBL = histMat;
binResRaw = Resolution*1e-12;      % binRes in seconds
binRes = downSampFactor*binResRaw;
fr = SyncRate;                  % Repetition frequency (in Hz)
tr = 1/fr;                      % Repetition period (in seconds)
numBins = floor(tr/binRes);
numBinsRaw = floor(tr/binResRaw);

histMatBL = zeros(numY, numX, numBins);
for kk = 1:downSampFactor
    dsTemp = downsample(transpose(reshape(histMat,[numY*numX,numBinsRaw])),downSampFactor,kk-1);
    histMatBL = histMatBL + reshape(dsTemp(1:numBins,:)',numY,numX,numBins);
end

numDetsBL = sum(histMat,3);
meanNumDetsBL = mean(numDetsBL(:));
medianNumDetsBL = median(numDetsBL(:));
maxNumDetsBL = max(numDetsBL(:));
disp(['Baseline dataset: Mean counts per pixel = ' num2str(meanNumDetsBL)]);
disp(['Baseline dataset: Median counts per pixel = ' num2str(medianNumDetsBL)]);
disp(['Baseline dataset: Max counts per pixel = ' num2str(maxNumDetsBL)]);

%% Load, preprocess pulse shape

fname1 = 'laser_pulse_shape_100_2018_12_06.txt';
pname1 = 'Data/';
fid = fopen([pname1 fname1]);
C1 = textscan(fid, '%f %f','HeaderLines', 8);
C2 = textscan(fid, '%f %f','HeaderLines', 1);
fclose(fid);

pulseShape = C2{1}(85:348);
binResP = C1{1}*1e-9;
numBinsP = length(pulseShape);

histResFactor = round(binRes/binResP);

% Pad pulse shape to correct length
pulseShapePad = [pulseShape; zeros(histResFactor*numBins-numBinsP,1)];

% Smooth pulse shape:
pulseShapeFilt = filtfilt(gausswin(3)/sum(gausswin(3)),1,pulseShapePad);

% Match data binning
sigPulse = zeros(numBins,1);
for ii = 1:histResFactor
    sigPulse = sigPulse+downsample(pulseShapeFilt,histResFactor,ii-1);
end
sigPulse = max(sigPulse-sigPulse(1),0);
sigPulse = sigPulse/sum(sigPulse);

%% Compute PDFs for Quantized S values
maxSigBL = max(SigBL(:));
minSigBL = min(SigBL(:));
rangeSigBL = maxSigBL-minSigBL;
numAlphas = 8;
alphaValsBL = min(SigBL(:))+(1:2:2*numAlphas-1)/(2*numAlphas)*rangeSigBL;

alphaBL_PDFs = zeros(numBins,numAlphas);
lmAlphaBL_PDFs = zeros(numBins,numAlphas);
for ii = 1:numAlphas
    alphaBL_PDFs(:,ii) = alphaValsBL(ii)*sigPulse+BgndBL/numBins;
    lmAlphaBL_PDFs(:,ii) = log(alphaBL_PDFs(:,ii));
end

%% Depth estimation for each pixel
Vxmin = -3; Vxmax = 3;
Vymin = -6; Vymax = 3;
Vxplot = linspace(Vxmin,Vymax,numX);
Vyplot = fliplr(linspace(Vymin,1.4,numY));
[Vxgrid,Vygrid] = meshgrid(Vxplot,Vyplot);

zBL = zeros(numY,numX);
alphaValIndex = zeros(numY,numX);
parfor ii = 1:numY
    for jj = 1:numX
        % get correct alphaVal index
        alphaValIndex(ii,jj) = min(1+floor(numAlphas*(SigBL(ii,jj)-minSigBL)/(rangeSigBL+eps)),numAlphas);
        xCor = cconv(squeeze(histMatBL(ii,jj,:))',flipud(lmAlphaBL_PDFs(:,alphaValIndex(ii,jj))),numBins);
        [~,maxLag] = max(xCor);
        zBL(ii,jj) = maxLag*binRes*1.5e8;
    end
end

%%

figPosZ = [300 300 350 400];
figPosL = [300 300 300 450];

rangeRefl = [0,6];
hfigBL_L = figure('position',figPosL); imagesc(SigBL,rangeRefl); axis image; 
colormap gray; colorbar('southoutside');
set(gca,'xtick',[],'ytick',[],'fontsize',18);
if printAllFigs
    fnameBL_L = [acqDate '_BL_L_' dtime];
    set(hfigBL_L,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosL(3), figPosL(4)])
    print(hfigBL_L,[prname fnameBL_L],'-dpdf','-r0');
end

rangeDepth1 = [0,4];
hfigBLz = figure('position',figPosZ); imagesc(zBL,rangeDepth1); axis image; % title('Baseline Depth Estimate');
set(gca,'xtick',[],'ytick',[],'fontsize',18); colorbar; colormap(jet);
if printAllFigs
    fnameBLz = [acqDate '_BLz_' dtime];
    set(hfigBLz,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosZ(3), figPosZ(4)])
    print(hfigBLz,[prname fnameBLz],'-dpdf','-r0');
end

rangeDepth2 = [0.56,0.9];
zBLmed = medfilt2d(zBL); % zBLmed = medfilt2(zBL,'symmetric');
speckleMask = abs(zBLmed-zBL)<0.03;
nearMask = (zBL>rangeDepth2(1) & zBL<rangeDepth2(2) & speckleMask);

%% 3D Point Cloud
alpha = pi/90;
Xbls = zBL(speckleMask(:)).*sin(2*alpha*Vxgrid(speckleMask(:)));
Ybls = zBL(speckleMask(:)).*cos(2*alpha*Vxgrid(speckleMask(:))).*sin(2*alpha*Vygrid(speckleMask(:)));
Zbls = zBL(speckleMask(:)).*cos(2*alpha*Vxgrid(speckleMask(:))).*cos(2*alpha*Vygrid(speckleMask(:)));

View3D = [20 11]; cRangeBL = rangeRefl; zlim3D = [-0.4 0.15]; ptSize = 6; 
fig3Dpos = [300 300 330 350]; %fig3Dpos = [300 300 450 450];
hfigBL = figure('position',fig3Dpos); scatter3(Xbls,Zbls,Ybls,ptSize,SigBL(speckleMask(:)),'filled');
% hfigBL.Renderer='Painters';
colormap gray; caxis(cRangeBL); 
colorbar('FontSize',18,'position',[0.86, 0.20, 0.05,0.7]);
set(gca,'color','k','gridcolor','w','FontSize',18,'ytick',[0.6,0.8]);
axis equal; ylim(rangeDepth2); zlim(zlim3D); set(gcf,'color','w');
view(View3D); title('Ground Truth');

if printAllFigs
    fnameBL = [acqDate '_BL_' dtime];
    set(hfigBL,'PaperPositionMode','Auto','InvertHardcopy','off','PaperUnits','points','PaperSize',[fig3Dpos(3), fig3Dpos(4)])
    print(hfigBL,[prname fnameBL],'-dpdf','-r0');
end

%% -----------------------------------------------------------------------
% ---------------------- High-Flux Measurements --------------------------
% ------------------------------------------------------------------------
% BgndHF = 0.01;
dwellTimeHF = 1000;

load([pname acqDate '_dt_data_acq' num2str(hfAcq) '_Histograms_' num2str(dwellTimeHF) '.mat']);
LambdaHF = LambdaMat;

% Assume background is 10% of total flux for every pixel
SDR = 0.1;
% BgndHF = SDR*LambdaHF;
SigHF = (1-SDR)*LambdaHF;%max(0.001,LambdaHF-BgndHF);

histMatHF = zeros(numY, numX, numBins);
for kk = 1:downSampFactor
    dsTemp = downsample(transpose(reshape(histMat,[numY*numX,numBinsRaw])),downSampFactor,kk-1);
    histMatHF = histMatHF + reshape(dsTemp(1:numBins,:)',numY,numX,numBins);
end

numDetsHF = sum(histMatHF,3);
meanNumDetsHF = mean(numDetsHF(:));
medianNumDetsHF = median(numDetsHF(:));
maxNumDetsHF = max(numDetsHF(:));
disp(['HF dataset: Mean counts per pixel = ' num2str(meanNumDetsHF)]);
disp(['HF dataset: Median counts per pixel = ' num2str(medianNumDetsHF)]);
disp(['HF dataset: Max counts per pixel = ' num2str(maxNumDetsHF)]);

%%
rangeReflHF = [0,1.5];
rangeReflMC = [0,4.5];
hfigHF_L = figure('position',figPosL); imagesc(LambdaHF,rangeReflMC); axis image; colormap gray; colorbar('southoutside');
% title(['$\Lambda^{\rm HF}$, numFrames = ' num2str(dwellTimeHF)],'interpreter','latex');
set(gca,'xtick',[],'ytick',[],'fontsize',18);
if printAllFigs
    fnameHF_L = [acqDate '_HF_L_' num2str(dwellTimeHF) '_' dtime];
    set(hfigHF_L,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosL(3), figPosL(4)])
    print(hfigHF_L,[prname fnameHF_L],'-dpdf','-r0');
end

countFracHF = sum(histMatHF,3)/dwellTimeHF;
hfigHF_num = figure('position',figPosL); imagesc(countFracHF,rangeReflHF); axis image; colormap gray; colorbar('southoutside');
% title('Num. Counts');
set(gca,'xtick',[],'ytick',[],'fontsize',18);
if printAllFigs
    fnameHF_num = [acqDate '_HF_num_' num2str(dwellTimeHF) '_' dtime];
    set(hfigHF_num,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosL(3), figPosL(4)])
    print(hfigHF_num,[prname fnameHF_num],'-dpdf','-r0');
end

figPosShort = [300 300 400 250];
hfigSat = figure('position',figPosShort); plot(LambdaHF(:),countFracHF(:),'o'); 
% hold on; plot(0:1,0:1,'linewidth',2)
xlabel('$\hat{\Lambda}^{\rm DML}$', 'interpreter', 'latex'); 
ylabel('$\hat{\Lambda}^{\rm AML}$', 'interpreter', 'latex'); 
set(gca,'fontsize',18);
% legend('Flux Estimate','Unit Slope','Location','best'); 
title('Flux Estimate Saturation');
if printAllFigs
    fnameHF_sat = [acqDate '_HF_sat' dtime];
    set(hfigSat,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosShort(3), figPosShort(4)])
    print(hfigSat,[prname fnameHF_sat],'-dpdf','-r0');
end

%% Compute PDFs for Quantized Lambda values
maxLamHF = max(LambdaHF(isfinite(LambdaHF)));
minLamHF = min(LambdaHF(isfinite(LambdaHF)));
rangeLamHF = maxLamHF-minLamHF;
alphaValsHF = min(SigHF(:))+(1:2:2*numAlphas-1)/(2*numAlphas)*rangeLamHF;

alphaHF_PDFs = zeros(numBins,numAlphas);
lmAlphaHF_PDFs = zeros(numBins,numAlphas);

Tbin  = binRes*1e9;
TdEst = td;
TeEst = te;
Tr = tr*1e9;

parfor ii = 1:numAlphas
    % Assume background is SDR fraction of total flux for every pixel
    BgndHF = SDR*alphaValsHF(ii);
    lam = (1-SDR)*alphaValsHF(ii)*sigPulse/Tbin+BgndHF/Tr;
    alphaHF_PDFs(:,ii) = fcn_PDF_PileupTdTe(Tbin, lam, Tr, TdEst, TeEst);
    lmAlphaHF_PDFs(:,ii) = log(alphaHF_PDFs(:,ii));
end

%% Depth estimation for each pixel
zMCPDF = zeros(numY,numX);
alphaValIndex = zeros(numY,numX);
tic;
parfor ii = 1:numY
    for jj = 1:numX
        % get correct alphaVal index
        alphaValIndex(ii,jj) = min(1+floor(numAlphas*(LambdaHF(ii,jj)-minLamHF)/(rangeLamHF+eps)),numAlphas);
        
        xCor = cconv(squeeze(histMatHF(ii,jj,:))',flipud(lmAlphaHF_PDFs(:,alphaValIndex(ii,jj))),numBins);
        [~,maxLag] = max(xCor);
        zMCPDF(ii,jj) = maxLag*binRes*1.5e8;
    end
end
toc;
%% Plots
errMCPDF = zMCPDF-zBL;
RMSE_MCPDF = sqrt(mean(errMCPDF(speckleMask).^2));
MAE_MCPDF = mean(abs(errMCPDF(speckleMask)));
RMSE_MCPDF_near = sqrt(mean(errMCPDF(nearMask).^2));
MAE_MCPDF_near = mean(abs(errMCPDF(nearMask)));

hfigMCPDFerr = figure('position',figPosZ); imagesc(errMCPDF,[-0.025 0.025]); axis image;
colormap jet; colorbar('Ticks',-0.02:0.01:0.02); set(gca,'FontSize',16,'YTick',[],'XTick',[]);
title(['MCPDF: MAE = ' num2str(MAE_MCPDF) ' m, RMSE = ' num2str(RMSE_MCPDF) ' m' ],'FontSize',10);

if printAllFigs
    fnameMCPDFerr = [acqDate '_MCPDFerr_' dtime];
    set(hfigMCPDFerr,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosZ(3), figPosZ(4)])
    print(hfigMCPDFerr,[prname fnameMCPDFerr],'-dpdf','-r0');
end

hfigMCPDFz = figure('position',figPosZ); imagesc(zMCPDF,rangeDepth1); axis image;
% title(['MCPDF Depth, Dwell Time ' num2str(dwellTimeHF) ' s, RMSE = ' num2str(RMSE_MCPDF) ' m']);
set(gca,'xtick',[],'ytick',[],'fontsize',18); colorbar; colormap(jet);
if printAllFigs
    fnameMCPDFz = [acqDate '_MCPDFz_' num2str(dwellTimeHF) '_' dtime];
    set(hfigMCPDFz,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosZ(3), figPosZ(4)])
    print(hfigMCPDFz,[prname fnameMCPDFz],'-dpdf','-r0');
end

% 3D Point Cloud
alpha = pi/90;
Xmcpdfs = zMCPDF(:).*sin(2*alpha*Vxgrid(:));
Ymcpdfs = zMCPDF(:).*cos(2*alpha*Vxgrid(:)).*sin(2*alpha*Vygrid(:));
Zmcpdfs = zMCPDF(:).*cos(2*alpha*Vxgrid(:)).*cos(2*alpha*Vygrid(:));

cRangeHF = rangeReflMC;
hfigMCPDF = figure('position',fig3Dpos); scatter3(Xmcpdfs,Zmcpdfs,Ymcpdfs,ptSize,LambdaHF(:),'filled');
% hfigMCPDF.Renderer='Painters';
colormap gray; axis equal;  caxis(cRangeHF);
set(gca,'color','k','gridcolor','w','FontSize',18,'ytick',[0.6,0.8]);
ylim([0.5 0.9]); zlim(zlim3D);  %colorbar('southoutside','FontSize',18); 
colorbar('FontSize',18,'position',[0.86, 0.20, 0.05,0.7]);
set(gcf,'color','w');
view(View3D); title('MCPDF');

if printAllFigs
    fnameMCPDF = [acqDate '_MCPDF_' dtime];
    set(hfigMCPDF,'PaperPositionMode','Auto','InvertHardcopy','off','PaperUnits','points','PaperSize',[fig3Dpos(3), fig3Dpos(4)])
    print(hfigMCPDF,[prname fnameMCPDF],'-dpdf','-r0');
end
%% HF (naive) Depth estimation for each pixel
zHF = zeros(numY,numX);
alphaValIndex = zeros(numY,numX);

parfor ii = 1:numY
    for jj = 1:numX
        % get correct alphaVal index
        alphaValIndex(ii,jj) = min(1+floor(numAlphas*(LambdaHF(ii,jj)-minLamHF)/(rangeLamHF+eps)),numAlphas);
        
        xCor = cconv(squeeze(histMatHF(ii,jj,:))',flipud(log(sigPulse+BgndBL/numBins)),numBins);
        [~,maxLag] = max(xCor);
        zHF(ii,jj) = maxLag*binRes*1.5e8;
    end
end

%% Plots
errHF = zHF-zBL;
RMSE_HF = sqrt(mean(errHF(speckleMask).^2));
MAE_HF = mean(abs(errHF(speckleMask)));
RMSE_HF_near = sqrt(mean(errHF(nearMask).^2));
MAE_HF_near = mean(abs(errHF(nearMask)));

hfigHFerr = figure('position',figPosZ); imagesc(errHF,[-0.025 0.025]); axis image;
colormap jet; colorbar('Ticks',-0.02:0.01:0.02); set(gca,'FontSize',16,'YTick',[],'XTick',[]);
title(['HF: MAE = ' num2str(MAE_HF) ' m, RMSE = ' num2str(RMSE_HF) ' m' ],'FontSize',10);

if printAllFigs
    fnameHFerr = [acqDate '_HFerr_' dtime];
    set(hfigHFerr,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosZ(3), figPosZ(4)])
    print(hfigHFerr,[prname fnameHFerr],'-dpdf','-r0');
end

hfigHFz = figure('position',figPosZ); imagesc(zHF,rangeDepth1); axis image;
% title(['HF Depth, Dwell Time ' num2str(dwellTimeHF) ' s, RMSE = ' num2str(RMSE_HF) ' m']);
set(gca,'xtick',[],'ytick',[],'fontsize',18); colorbar; colormap(jet)
if printAllFigs
    fnameHFz = [acqDate '_HFz_' num2str(dwellTimeHF) '_' dtime];
    set(hfigHFz,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosZ(3), figPosZ(4)])
    print(hfigHFz,[prname fnameHFz],'-dpdf','-r0');
end


% 3D Point Cloud
alpha = pi/90;
Xhfs = zHF(:).*sin(2*alpha*Vxgrid(:));
Yhfs = zHF(:).*cos(2*alpha*Vxgrid(:)).*sin(2*alpha*Vygrid(:));
Zhfs = zHF(:).*cos(2*alpha*Vxgrid(:)).*cos(2*alpha*Vygrid(:));

hfigHF = figure('position',fig3Dpos); scatter3(Xhfs,Zhfs,Yhfs,ptSize,countFracHF(:),'filled');
% hfigHF.Renderer='Painters';
colormap gray; axis equal;  caxis(rangeReflHF);
set(gca,'color','k','gridcolor','w','FontSize',18,'ytick',[0.6,0.8]);
ylim([0.5 0.9]); zlim(zlim3D);  %colorbar('southoutside','FontSize',18); 
colorbar('FontSize',18,'position',[0.86, 0.20, 0.05,0.7]);
set(gcf,'color','w');
view(View3D); title('HF');
if printAllFigs
    fnameHF = [acqDate '_HF_' dtime];
    set(hfigHF,'PaperPositionMode','Auto','InvertHardcopy','off','PaperUnits','points','PaperSize',[fig3Dpos(3), fig3Dpos(4)])
    print(hfigHF,[prname fnameHF],'-dpdf','-r0');
end

%%
% ------------------------------------------------------------------------
% ---------------------- Low-Flux Measurements ---------------------------
% ------------------------------------------------------------------------

BgndLF = BgndBL;%1;
dwellTimeLF = dwellTimeHF; %*100; %1000;

load([pname acqDate '_dt_data_acq' num2str(lfAcq) '_Histograms_' num2str(dwellTimeLF) '.mat']);
LambdaLF = LambdaMat;
OD_LF = OD;
SigLF = max(0.001,LambdaLF-BgndLF);
figure('position',figPosL); imagesc(LambdaLF,[0,20]); axis image; colormap gray; colorbar;
title(['$\Lambda^{\rm LF}$, numFrames = ' num2str(dwellTimeLF)],'interpreter','latex');

% histMatLF = histMat;
histMatLF = zeros(numY, numX, numBins);
for kk = 1:downSampFactor
    dsTemp = downsample(transpose(reshape(histMat,[numY*numX,numBinsRaw])),downSampFactor,kk-1);
    histMatLF = histMatLF + reshape(dsTemp(1:numBins,:)',numY,numX,numBins);
end

numDetsLF = sum(histMatLF,3);
meanNumDetsLF = mean(numDetsLF(:));
medianNumDetsLF = median(numDetsLF(:));
maxNumDetsLF = max(numDetsLF(:));
disp(['LF dataset: Mean counts per pixel = ' num2str(meanNumDetsLF)]);
disp(['LF dataset: Median counts per pixel = ' num2str(medianNumDetsLF)]);
disp(['LF dataset: Max counts per pixel = ' num2str(maxNumDetsLF)]);

%%
countFracLF = sum(histMatLF,3)/dwellTimeLF;
fluxLF = countFracLF*10^OD_LF;
hfigLF_num = figure('position',figPosL); imagesc(fluxLF,rangeRefl); axis image; colormap gray; colorbar('southoutside');
% title('Num. Counts');
set(gca,'xtick',[],'ytick',[],'fontsize',18);
if printAllFigs
    fnameLF_num = [acqDate '_LF_num_' num2str(dwellTimeLF) '_' dtime];
    set(hfigLF_num,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosL(3), figPosL(4)])
    print(hfigLF_num,[prname fnameLF_num],'-dpdf','-r0');
end

maxLamLF = max(LambdaLF(isfinite(LambdaLF)));

if maxLamLF > 2*maxSigBL    % Outlier maximum estimates, use baseline values
    rangeSigLF = maxSigBL-minSigBL;
    alphaValsLF = min(SigLF(:))+(1:2:2*numAlphas-1)/(2*numAlphas)*rangeSigLF;
else    % Use data-driven values
    minLamLF = min(LambdaLF(isfinite(LambdaLF)));
    rangeLamLF = maxLamLF-minLamLF;
    alphaValsLF = min(SigLF(:))+(1:2:2*numAlphas-1)/(2*numAlphas)*rangeLamLF;
end

alphaLF_PDFs = zeros(numBins,numAlphas);
lmAlphaLF_PDFs = zeros(numBins,numAlphas);
for ii = 1:numAlphas
    alphaLF_PDFs(:,ii) = alphaValsLF(ii)*sigPulse+BgndLF/numBins;
    lmAlphaLF_PDFs(:,ii) = log(alphaLF_PDFs(:,ii));
end

clear histMat;

% Depth estimation for each pixel
zLF = zeros(numY,numX);
alphaValIndexLF = zeros(numY,numX);

parfor ii = 1:numY
    for jj = 1:numX
        % get correct alphaVal index
        alphaValIndexLF(ii,jj) = max(1,min(1+floor(numAlphas*(fluxLF(ii,jj)-minSigBL)/(rangeSigLF+eps)),numAlphas));
        
        xCor = cconv(squeeze(histMatLF(ii,jj,:))',flipud(lmAlphaLF_PDFs(:,alphaValIndexLF(ii,jj))),numBins);
        [~,maxLag] = max(xCor);
        zLF(ii,jj) = maxLag*binRes*1.5e8;
    end
end

% Plots

errLF = zLF-zBL;
RMSE_LF = sqrt(mean(errLF(speckleMask).^2));
MAE_LF = mean(abs(errLF(speckleMask)));
RMSE_LF_near = sqrt(mean(errLF(nearMask).^2));
MAE_LF_near = mean(abs(errLF(nearMask)));

hfigLFerr = figure('position',figPosZ); imagesc(errLF,[-0.025 0.025]); axis image;
colormap jet; colorbar('Ticks',-0.02:0.01:0.02); set(gca,'FontSize',16,'YTick',[],'XTick',[]);
title(['LF: MAE = ' num2str(MAE_LF) ' m, RMSE = ' num2str(RMSE_LF) ' m' ],'FontSize',10);
if printAllFigs
    fnameLFerr = [acqDate '_LFerr_' num2str(dwellTimeLF) '_' dtime];
    set(hfigLFerr,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosZ(3), figPosZ(4)])
    print(hfigLFerr,[prname fnameLFerr],'-dpdf','-r0');
end

hfigLFz = figure('position',figPosZ); imagesc(zLF,rangeDepth1); axis image;
title(['LF Depth, Dwell Time ' num2str(dwellTimeLF) ' s, RMSE = ' num2str(RMSE_LF) ' m']);
set(gca,'xtick',[],'ytick',[],'fontsize',18); colorbar; colormap(jet)
if printAllFigs
    fnameLFz = [acqDate '_LFz_' num2str(dwellTimeLF) '_' dtime];
    set(hfigLFz,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosZ(3), figPosZ(4)])
    print(hfigLFz,[prname fnameLFz],'-dpdf','-r0');
end


% 3D Point Cloud
alpha = pi/90;
Xhfs = zLF(:).*sin(2*alpha*Vxgrid(:));
Yhfs = zLF(:).*cos(2*alpha*Vxgrid(:)).*sin(2*alpha*Vygrid(:));
Zhfs = zLF(:).*cos(2*alpha*Vxgrid(:)).*cos(2*alpha*Vygrid(:));

hfigLF = figure('position',fig3Dpos); scatter3(Xhfs,Zhfs,Yhfs,ptSize,fluxLF(:),'filled');
% hfigLF.Renderer='Painters';
colormap gray; axis equal;  caxis(cRangeBL);
set(gca,'color','k','gridcolor','w','FontSize',18,'ytick',[0.6,0.8]);
ylim([0.5 0.9]); zlim(zlim3D);  %colorbar('southoutside','FontSize',18); 
colorbar('FontSize',18,'position',[0.86, 0.20, 0.05,0.7]);
set(gcf,'color','w');
view(View3D); title('LF');

if printAllFigs
    fnameLF = [acqDate '_LF_' num2str(dwellTimeLF) '_' dtime];
    set(hfigLF,'PaperPositionMode','Auto','InvertHardcopy','off','PaperUnits','points','PaperSize',[fig3Dpos(3), fig3Dpos(4)])
    print(hfigLF,[prname fnameLF],'-dpdf','-r0');
end

%% Load pre-processed MCHC results
load([pname 'MCHC_result.mat']);

errMCHC = zMCHC-zBL;
RMSE_MCHC = sqrt(mean(errMCHC(speckleMask).^2));
MAE_MCHC = mean(abs(errMCHC(speckleMask)));
RMSE_MCHC_near = sqrt(mean(errMCHC(nearMask).^2));
MAE_MCHC_near = mean(abs(errMCHC(nearMask)));

hfigMCHCerr = figure('position',figPosZ); imagesc(errMCHC,[-0.025 0.025]); axis image;
colormap jet; colorbar('Ticks',-0.02:0.01:0.02); set(gca,'FontSize',16,'YTick',[],'XTick',[]);
title(['MCHC: MAE = ' num2str(MAE_MCHC) ' m, RMSE = ' num2str(RMSE_MCHC) ' m' ],'FontSize',10);

if printAllFigs
    fnameMCHCerr = [acqDate '_MCHCerr_' dtime];
    set(hfigMCHCerr,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosZ(3), figPosZ(4)])
    print(hfigMCHCerr,[prname fnameMCHCerr],'-dpdf','-r0');
end

hfigMCHCz = figure('position',figPosZ); imagesc(zMCHC,rangeDepth1); axis image;
% title(['MCHC Depth, Dwell Time ' num2str(dwellTimeHF) ' s, RMSE = ' num2str(RMSE_MCHC) ' m']);
set(gca,'xtick',[],'ytick',[],'fontsize',18); colorbar; colormap(jet);
if printAllFigs
    fnameMCHCz = [acqDate '_MCHCz_' num2str(dwellTimeHF) '_' dtime];
    set(hfigMCHCz,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[figPosZ(3), figPosZ(4)])
    print(hfigMCHCz,[prname fnameMCHCz],'-dpdf','-r0');
end

% 3D Point Cloud
alpha = pi/90;
XMCHCs = zMCHC(:).*sin(2*alpha*Vxgrid(:));
YMCHCs = zMCHC(:).*cos(2*alpha*Vxgrid(:)).*sin(2*alpha*Vygrid(:));
ZMCHCs = zMCHC(:).*cos(2*alpha*Vxgrid(:)).*cos(2*alpha*Vygrid(:));

cRangeHF = [0,4.5];
hfigMCHC = figure('position',fig3Dpos); scatter3(XMCHCs,ZMCHCs,YMCHCs,ptSize,LambdaHF(:),'filled');
colormap gray; axis equal;  caxis(cRangeHF);
set(gca,'color','k','gridcolor','w','FontSize',18,'ytick',[0.6,0.8]);
ylim([0.5 0.9]); zlim(zlim3D);  %colorbar('southoutside','FontSize',18); 
colorbar('FontSize',18,'position',[0.86, 0.20, 0.05,0.7]);
set(gcf,'color','w');
view(View3D); title('MCHC');

if printAllFigs
    fnameMCHC = [acqDate '_MCHC_' dtime];
    set(hfigMCHC,'PaperPositionMode','Auto','InvertHardcopy','off','PaperUnits','points','PaperSize',[fig3Dpos(3), fig3Dpos(4)])
    print(hfigMCHC,[prname fnameMCHC],'-dpdf','-r0');
end

