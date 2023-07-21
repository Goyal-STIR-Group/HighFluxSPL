function [piHat,trans] =  fcn_PDF_Pileup0_pts2pts_Te(bins,Tbin, lam, Tr, Td, Te, flag_piHat)

numBins = length(bins);
piHat = zeros(numBins,1);

nTd = Td/Tbin;
nTe = Te/Tbin;
nTr = Tr/Tbin;

int_lam = cumsum(lam)*Tbin;
Lam = sum(lam)*Tbin;
lam_ext = [lam,lam,lam];
int_lam_ext = cumsum(lam_ext)*Tbin;

hx = zeros(size(bins));
for i = 1:numBins
    exponent = sum(lam_ext(i+1:i+nTd))*Tbin;
    hx(i) = lam(i)*exp(exponent);
end
hx_ext = [hx,hx,hx];

int_hx = zeros(size(hx));
for i = 1:numBins
    int_hx(i) = sum(hx_ext(i+nTd+1:i+nTe))*Tbin;  
end

trans = zeros(numBins);
for i = 1:numBins
    if i  <= nTr - nTe
        exponent = int_lam + Lam -int_lam(i+nTd) ;
        exponent(i + nTe + 1: end) = exponent(i + nTe + 1: end)  - Lam;
        lam_exp =  lam.*exp(-exponent);
                
        ind1 = [1: i + 2*nTd - nTr, i + nTe+1: nTr];
        trans(i,ind1) =lam_exp(ind1)*(1+exp(-Lam)*int_hx(i))/(1- exp(-Lam));
        
        ind2 = i + 2*nTd - nTr +1: i + nTd + nTe - nTr;
        for j = 1:length(ind2)
            trans(i,ind2(j)) =lam_exp(ind2(j))*( (1+exp(-Lam)*int_hx(i))/(1- exp(-Lam)) + sum(hx(i+nTd+1: ind2(j)+nTr-nTd))*Tbin);
        end
        
        ind3 =  i + nTd + nTe - nTr + 1: i + nTe;
        trans(i,ind3) =lam_exp(ind3)*(1+int_hx(i))/(1- exp(-Lam));
        
    elseif i <= 2*nTr - nTd - nTe
        
        exponent = int_lam_ext(2*nTr + 1: 2*nTr + numBins) - int_lam_ext(i + nTd);
        exponent(i + nTe - nTr + 1:end) =  exponent(i + nTe - nTr +1 :end) - Lam;
        lam_exp =  lam.*exp(-exponent);
        
        ind1 = [1:i+nTe-nTr, i + nTd + nTe - nTr+1: nTr];
        trans(i,ind1) = lam_exp(ind1)* (1+int_hx(i))/(1-exp(-Lam));
        
        ind2 = i + nTe - nTr + 1: i + 2*nTd - nTr;
        trans(i,ind2) = lam_exp(ind2)* (1+exp(-Lam)*int_hx(i))/(1-exp(-Lam));
        
        ind3 = i + 2*nTd - nTr+1: i+nTd+nTe-nTr;
        for j = 1:length(ind3)
            trans(i,ind3(j)) = lam_exp(ind3(j))* ((1+exp(-Lam)*int_hx(i))/(1-exp(-Lam)) + sum(hx_ext(i+nTd+1:ind3(j)+nTr-nTd))*Tbin);
        end
        
    elseif i <= 2*nTr - 2*nTd
        exponent = int_lam_ext(2*nTr + 1: 2*nTr + numBins) - int_lam_ext(i + nTd);
        exponent(i + nTe - nTr + 1:end) =  exponent(i + nTe - nTr + 1:end) - Lam;
        lam_exp =  lam.*exp(-exponent);
        
        ind1 = 1:i+nTd+nTe-2*nTr;
        for j = 1:length(ind1)
            trans(i,ind1(j)) = lam_exp(ind1(j))* ((1+exp(-Lam)*int_hx(i))/(1-exp(-Lam)) + sum(hx_ext(i+nTd+1:ind1(j) + 2*nTr - nTd))*Tbin);
        end
        
        ind2 = i + nTe +nTd - 2*nTr + 1: i + nTe - nTr;
        trans(i,ind2) = lam_exp(ind2)* (1+int_hx(i))/(1-exp(-Lam));
        
        ind3 = i + nTe - nTr + 1: i+2*nTd-nTr;
        trans(i,ind3) = lam_exp(ind3)* (1+exp(-Lam)*int_hx(i))/(1-exp(-Lam));
        
        ind4 = i+2*nTd-nTr + 1: nTr;
        for j = 1:length(ind4)
            trans(i,ind4(j)) = lam_exp(ind4(j))* ((1+exp(-Lam)*int_hx(i))/(1-exp(-Lam)) + sum(hx_ext(i+nTd+1: ind4(j)+nTr-nTd))*Tbin);
        end
        
    else
        
        exponent = int_lam_ext(2*nTr + 1: 2*nTr + numBins) - int_lam_ext(i + nTd);
        exponent(i + nTe - nTr + 1:end) =  exponent(i + nTe - nTr + 1:end) - Lam;
        lam_exp =  lam.*exp(-exponent);
        
        ind1 = [1:i+2*nTd-2*nTr, i+nTe-nTr+1:nTr];
        trans(i,ind1) = lam_exp(ind1)* (1+exp(-Lam)*int_hx(i))/(1-exp(-Lam));
        
        ind2 = i+2*nTd-2*nTr+1: i + nTe + nTd-2*nTr;
        for j = 1:length(ind2)
            trans(i,ind2(j)) = lam_exp(ind2(j))* ((1+exp(-Lam)*int_hx(i))/(1-exp(-Lam)) + sum(hx_ext(i+nTd+1:ind2(j) + 2*nTr - nTd))*Tbin);
        end
        
        ind3 = i + nTe + nTd-2*nTr+1:i+nTe-nTr;
        trans(i,ind3) = lam_exp(ind3)* (1+int_hx(i))/(1-exp(-Lam));
    end
     trans(i,:) = trans(i,:)/sum(trans(i,:));
end

%% Compute stationary density function
if flag_piHat
    [f,~] = eigs(trans',1); 

    if f(1)<0
        f = -f;
    end
    if any(f<0)
        error('Transition matrix is wrong.')
    end
    piHat = f/sum(f);
end
