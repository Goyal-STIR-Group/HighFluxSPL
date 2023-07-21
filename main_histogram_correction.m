%% Simulating electronics and detector dead time
% Joshua Rapp
% Boston University
%
% Markov Chain-based Histogram Correction accounting for both detector and 
% electronics dead times (Td and Te).
%
% Generates Fig. 6 in the paper.
%
% Please cite
% J. Rapp, Y. Ma, R. Dawson, and V. Goyal, 
% “High-flux single-photon lidar,” Optica, vol. 8, no. 1, 
% pp. 30–39, Nov. 2021, doi: 10.1364/optica.403190.
%
% Also compares against
% S. Isbaner et al., “Dead-time correction of fluorescence lifetime 
% measurements and fluorescence lifetime imaging,” Optics Express, 
% vol. 24, no. 9, p. 9429, 2016, doi: 10.1364/OE.24.009429.

clear; close all; clc;
addpath(genpath('utils'));

% addpath(genpath('DeadtimeCorrectionCode'))

% Signal parameters
sigP = 2;       % pulse std [ns]
Tr = 100;       % illumination period [ns]
Td =55;         % detector dead time [ns]
Te = 80;        % electronics dead time [ns]
Ztrue =50;      % depth [ns]
Tbin = 0.05;    % bin size [ns]
bins = Tbin/2 : Tbin : Tr-Tbin/2; % bin centers [ns]
numBins = length(bins);

% Set mean Signal and Background levels
S_all = 2; %[1,2];
B_all = 2; %[1,2];

groundtruth = zeros(numBins,length(S_all),length(B_all));
measurement = zeros(numBins,length(S_all),length(B_all));
result_Proposed = zeros(numBins,length(S_all),length(B_all));
result_Isbaner = zeros(numBins,length(S_all),length(B_all));

for iS = 1:length(S_all)
    for iB = 1:length(B_all)       
        S = S_all(iS);
        B = B_all(iB);
        Lam = S+B;
        
        disp(['Processing S=' num2str(S) ', B=' num2str(B) '...']);
        
        % Ground truth photon arrival distribution
        lam_true =S*normpdf(bins,Ztrue,sigP)+B/Tr;
        groundtruth(:,iS,iB) = lam_true(:);
        f_arriv_true = lam_true/sum(lam_true*Tbin);
        flag_piHat = 1;
        
        % Detection distribution affected by dead times
        f_true =  fcn_PDF_Pileup0_pts2pts_Te(bins,Tbin, lam_true, Tr, Td, Te,flag_piHat);
        f_true = f_true/sum(f_true*Tbin);
        measurement(:,iS,iB) = f_true;
        
        % solve the inverse problem via proposed MCHC method
        max_iter = 100;
        lam_est = MarkovChainHistCorrection_Te(f_true,Lam,bins,Tbin,Tr,Td,Te,max_iter);
        result_Proposed(:,iS,iB) = lam_est;
        
        % solve the inverse problem using Isbaner's iterative idea
        lam_est_Isbaner = DeadTimeCorrection(f_true,Td,Te,Lam);
        result_Isbaner(:,iS,iB) = lam_est_Isbaner/sum(lam_est_Isbaner*Tbin)*Lam;
    end
end

%% Plots
cmap = lines(5);
figure('color','white','units','centimeters','position',[2 2 16 14]);
for iS = 1:length(S_all)
    for iB = 1:length(B_all)
        S = S_all(iS);
        B = B_all(iB);
        subplot(length(S_all)+1,length(B_all),(iS-1)*length(B_all)+iB)
        semilogy(bins,groundtruth(:,iS,iB),'black-','linewidth',1.5);
        hold on;semilogy(bins,measurement(:,iS,iB)/sum(measurement(:,iS,iB)*Tbin)*Lam,'-','linewidth',1.5,'color',cmap(1,:));
        hold on;semilogy(bins,result_Proposed(:,iS,iB),'--','linewidth',1.5,'color',cmap(2,:));
        hold on;semilogy(bins,result_Isbaner(:,iS,iB),'--','linewidth',1.5,'color',cmap(5,:));
        hold off;
        title(['S = ',num2str(S),', B = ',num2str(B)]);
        xlabel('Time [ns]'); ylim([1e-3,1]);
        if iS == length(S_all) && iB == length(B_all)
            legend('Ground truth','Measurement (scaled)','Proposed','Isbaner',...
                'location','southoutside','orientation','horizontal')
        end
    end
end
