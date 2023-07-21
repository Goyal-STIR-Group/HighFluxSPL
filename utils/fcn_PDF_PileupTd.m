function [ piHat ] = fcn_PDF_PileupTd(lam,Td,Tr,Tbin)
%FCN_PDF_PILEUPTD predicts the distortion to the photon detection
%distribution with 
% - arrival intensity lam, 
% - repetition period Tr
% - detector dead time Td
% - time bin Tbin
%
% Code written by Yanting Ma
%
% Please cite J. Rapp, Y. Ma, R. M. A. Dawson, and V. K. Goyal, 
% “Dead time compensation for high-flux ranging,” 
% IEEE Transactions on Signal Processing, pp. 1–1, Oct. 2019, 
% doi: 10.1109/TSP.2019.2914891.


nTd = floor(mod(Td,Tr)/Tbin);
nTr = floor(Tr/Tbin);
numBin = length(lam);

int_lam = cumsum(lam)*Tbin;
Lam = sum(lam)*Tbin;

trans = zeros(numBin);
for i = 1:numBin
    if i + nTd <= nTr
        exponent = int_lam + Lam -int_lam(i+nTd) ;
        exponent(i + nTd + 1: end) = exponent(i + nTd + 1: end)  - Lam;
        trans(i,:) = lam.*exp(-exponent)/(1- exp(-Lam));
        trans(i,:) = trans(i,:)/sum(trans(i,:));
    else
       exponent = int_lam + Lam -int_lam(i+nTd - nTr) ;
       exponent(i + nTd -nTr + 1: end) = exponent(i + nTd - nTr + 1: end)  - Lam;
       trans(i,:) = lam.*exp(-exponent)/(1- exp(-Lam));
       trans(i,:) = trans(i,:)/sum(trans(i,:));
    end
end

[f,~] = eigs(trans',1); 

if f(1)<0
    f = -f;
end
if any(f<0)
    error('Transition matrix is wrong.')
end
piHat = f/sum(f*Tbin);