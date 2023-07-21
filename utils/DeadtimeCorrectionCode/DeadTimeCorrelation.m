function err = DeadTimeCorrelation(p,t,y)

% function for calculating the correlation function for a cw-source in a
% TCSPC measurement setup
%
% input parameters:
% p - parameter vector: [electronics deadtime, detector deadtime, photon hit rate]
% (photon hit rate is average number of photons per time unit hitting the
% detector with supposed quantum yield of detection equal to 1)
% t - correlation time axis
% y (optional) - measured correlation curve: if this input is given, the
% function returns the fit error, otherwise the computed correlation curve
% itself

t = t(:);
deadtime = round(p(1)+1);
apd = round(p(2));
eps = p(3);

if apd==0
    decay0 = exp(-eps*(t-0.5));
    decay0 = decay0/sum(decay0);
    decay = decay0;
    jmax = floor(t(end)/deadtime);
    z = zeros(numel(t),1);
    for j=1:jmax
        z(j*deadtime:end) = z(j*deadtime:end) + decay(1:end-j*deadtime+1);
        decay = conv(decay,decay0);
        decay = decay(1:numel(t));
    end
    z = z/eps;
else
    tmp = DeadTimeCorrelation([apd 0 eps],t);
    decay0 = tmp(deadtime+apd)*exp(-eps*(t-apd+0.5));
    decay0(t<=apd) = tmp(t<=apd+deadtime & t>deadtime);
    decay0 = decay0/sum(decay0);
    decay = decay0;
    jmax = floor(t(end)/deadtime);
    z = zeros(numel(t),1);
    for j=1:jmax
        z(j*deadtime:end) = z(j*deadtime:end) + decay(1:numel(t)-j*deadtime+1);
        decay = conv(decay,decay0);
        decay = decay(1:numel(decay0));
    end
    z = z/eps;
end


if nargin<3
    err = z;
else
    c = z\y(:);
    z = c*z;
    plot(t,y,'o',t,z); drawnow
    err = sum((y-z).^2)
end

return

plot(autotime,tmp)
patch([0 2^6 2^6 0],[0 0 1 1],'b','facealpha',0.1)
patch([2^6 2^6+2^4 2^6+2^4 2^6],[0 0 1 1],'y','facealpha',0.1)
xlabel('correlation time'); ylabel('autocorrelation');
for j=1:6 s{j}=['\epsilon\itP\rm = ' mnum2str(epsv(j),1,1)]; end
legend(s,4)