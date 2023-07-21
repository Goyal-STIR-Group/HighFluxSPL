function [ LFDetTimes,HFDetTimes ] = ...
    fcn_gen_DE_Dets3( Nr,S,B,Tr,Td,Ztrue,sigP,Te,Te_subMean,Td_sig,detDecay,LFflag)
%FCN_GEN_DE_DETS3 generates a sequence of photon arrivals and then produces
% (1) A Bernoulli-thinned sequence meeting the desired low-flux condition,
% and
% (2) A deterministically-culled sequence accounting for dead time effects
%
% Nr = base number illuminations
% S = expected number of signal arrivals per period
% B = expected number of background arrivals per period
% Tr = illumination period (assumed nanoseconds)
% Td = detector dead time (assumed ns)
% Ztrue = true signal mean
% sigP = standard deviation of signal pulse
% Te = electronics dead time
% Te_subMean = mean of exponential RV subtracted from Te, must be
%       nonnegative
% Td_sig = standard deviation of detector dead time, must be nonnegative
% detDecay = decay of exponential detector response function
% LFflag = produce low-flux sequence? 1 if yes, 0 if no.
%
% Code written by Joshua Rapp
%
% Please cite J. Rapp, Y. Ma, R. Dawson, and V. Goyal, 
% “High-flux single-photon lidar,” Optica, vol. 8, no. 1, 
% pp. 30–39, Nov. 2021, doi: 10.1364/optica.403190.

numEnough = Nr*0.05;

warm_up = 1000; % number of illuminations starting from the beginning to be discarded;
Nr = ceil(1.1*Nr) + warm_up;

preCheck = 2*ceil((S+B)*(Td+Te)/Tr);

enoughDets = 0;
while ~enoughDets
    %% Generate Photon Arrivals
    numSig = poissrnd(S,Nr,1);   % Number of signal arrivals in period ii
    numBack = poissrnd(B,Nr,1);  % Number of background arrivals in period ii
    totNumSig = sum(numSig);
    totNumBack = sum(numBack);
    numArriv = totNumSig+totNumBack;
    arrivTimes = zeros(numArriv,1);
    
    rawSigTimes = Ztrue+sigP*randn(totNumSig,1);    % Gaussian sig photon arrival times
    rawBackTimes = Tr*rand(totNumBack,1);           % background photon arrival times
    cumNumSig = [0; cumsum(numSig)];
    cumNumBack = [0; cumsum(numBack)];
    numDets = numSig+numBack;
    
    indexNr = 1;
    for ii = 1:Nr
        arrivTimes(indexNr:indexNr+numDets(ii)-1) = ...
            sort([rawSigTimes(cumNumSig(ii)+1:cumNumSig(ii+1)); ...
            rawBackTimes(cumNumBack(ii)+1:cumNumBack(ii+1))])+(ii-1)*Tr;
        indexNr = indexNr+numDets(ii);
    end
    
    if LFflag == 1
        %% Low-flux photon detections
        % First attenuate overall signal so that arrival rate is approximately 5%
        % of illumination rate. Equivalent to placing neutral density filter in
        % front of detector.
        
        PrDet = min(0.05/(S+B),1);
        numLFArriv = binornd(numArriv,PrDet);
        LFarrivTimes = sort(randsample(arrivTimes,numLFArriv));
        
        LFDetTimes = dtDets(LFarrivTimes,Td,Te,Te_subMean,Td_sig,detDecay,preCheck);
        
        % Discard "warm-up" detections (before stationary distribution)
        LFDetTimes = LFDetTimes(LFDetTimes>warm_up*Tr) - warm_up*Tr;
        
        % Enough low-flux detections for comparison?
        numLF = length(LFDetTimes);
        if numLF > numEnough
            enoughDets = 1;
        else
            numNeeded = numEnough-numLF;
            Nr = Nr + numNeeded/0.01;
        end
    else
        enoughDets = 1;
        LFDetTimes = [];
    end
    
end

%% Generate Dead-time-affected Photon Detections
HFDetTimes = dtDets(arrivTimes,Td,Te,Te_subMean,Td_sig,detDecay,preCheck);

% Discard "warm-up" detections (before stationary distribution)
HFDetTimes = HFDetTimes(HFDetTimes>warm_up*Tr) - warm_up*Tr;

end

function dtDetTimes = dtDets(arrivTimes,Td,Te,Te_subMean,Td_sig,detDecay,preCheck)
numArriv = length(arrivTimes);

dtDetTimes = zeros(numArriv,1);
TeActual = Te-exprnd(Te_subMean,numArriv,1);
TdActual = Td+Td_sig*randn(numArriv,1);
dlayArrivTimes = arrivTimes+exprnd(detDecay,numArriv,1);


jj = 1; % arrival index
kk = 0; % detection index
while jj <= numArriv
    kk = kk + 1;
    dtDetTimes(kk) = dlayArrivTimes(jj);
    
    lastDIndx = jj; % Last detector event is detected photon
    
    % In order to find next detected arrival time, search increasing larger
    % sets of subsequent arrival times (starts small to speed up algorithm)
    nextCheck = preCheck;
    findFlag = 1;
    while nextCheck < numArriv && findFlag
        checkMax = min(lastDIndx+nextCheck,numArriv);
        checkTimes = dlayArrivTimes(lastDIndx+1:checkMax);
        nextDIndx = find(checkTimes>dlayArrivTimes(lastDIndx)+TdActual(lastDIndx),1);
        if isempty(nextDIndx)
            nextCheck = nextCheck*10;
        else
            findFlag = 0;
        end
    end
    
    % If detector event occurred during electronics dead time:
    while dlayArrivTimes(lastDIndx + nextDIndx)-dtDetTimes(kk)< TeActual(jj)
        lastDIndx = lastDIndx + nextDIndx;
        
        % Find next detector event since last detector event:
        %         nextDIndx = find(dlayArrivTimes(lastDIndx+1:end)>dlayArrivTimes(lastDIndx)+TdActual(lastDIndx),1);
        nextCheck = preCheck;
        findFlag = 1;
        while nextCheck < numArriv && findFlag
            checkMax = min(lastDIndx+nextCheck,numArriv);
            checkTimes = dlayArrivTimes(lastDIndx+1:checkMax);
            nextDIndx = find(checkTimes>dlayArrivTimes(lastDIndx)+TdActual(lastDIndx),1);
            if isempty(nextDIndx)
                nextCheck = nextCheck*10;
            else
                findFlag = 0;
            end
        end
    end
    
    % When both detector and electronics are reset:
    jj = lastDIndx + nextDIndx;
end
dtDetTimes = dtDetTimes(1:kk);
end