function lam_est = MarkovChainHistCorrection_Te(f_est,Lam,bins,Tbin,Tr,Td,Te,max_iter)
%MarkovChainHistCorrection_Te performs MCHC method
% - distorted histogram f_est
% - arrival rate Lam, 
% - time bins
% - repetition period Tr
% - detector dead time Td
% - electronics dead time Te
% - time bin Tbin
% where Td < Te <= Tr <= 2Td.
%
% Code written by Yanting Ma
%
% Please cite J. Rapp, Y. Ma, R. Dawson, and V. Goyal, 
% “High-flux single-photon lidar,” Optica, vol. 8, no. 1, 
% pp. 30–39, Nov. 2021, doi: 10.1364/optica.403190.

f_est = f_est(:)';
f_est = f_est/sum(f_est*Tbin); 
lam_est = f_est*Lam; % initialization


flag_piHat = 0;
for iter = 1:max_iter    
    [~,P] =  fcn_PDF_Pileup0_pts2pts_Te(bins,Tbin, lam_est, Tr, Td, Te,flag_piHat);
    s_est = (f_est*P)./(lam_est+realmin); 
    lam_est_prev = lam_est;
    lam_est = f_est./(s_est+ realmin);
    lam_est = lam_est/sum(lam_est*Tbin)*Lam;
    
    invErr(iter) = norm(lam_est_prev - lam_est)/norm(lam_est_prev);
    if norm(lam_est_prev - lam_est)/norm(lam_est_prev) < 1e-5
        break
    end
end

if norm(lam_est_prev - lam_est)/norm(lam_est_prev) > 1e-1
    %lam_est = f_est*Lam;
    lam_est = lam_est(:);
    warning('MCHC failed.')
else
    lam_est = lam_est(:);
end

figure; plot(invErr);