function [z, xv] = mHist(x,xv,weight)

if size(x,1)==1 || size(x,2)==1
    x = x(:);
    if nargin>2
        weight = weight(:);
        weight(~isfinite(x)) = [];
    end
    x(~isfinite(x)) = [];
end
if nargin>1 && ~isempty(xv)
    if size(xv)==1
        n = xv;
    else
        n = [];
        xmin = xv(1);
        xmax = xv(end);
        ind = x>xmax | x<xmin;
        x(ind) = [];
        if nargin>2
            weight(ind) = [];
        end
        if length(unique(diff(xv)))==1
            dx = xv(2)-xv(1);
            x = round((x-xmin)/dx)+1;
            xmax = round((xmax-xmin)/dx)+1;
        else
            x = round(interp1(xv,1:length(xv),x));
            xmax = round(interp1(xv,1:length(xv),xmax));
        end
    end
else
    n = 100;
end
if ~isempty(n)
    xmin = min(min(x));
    xmax = max(max(x));
    dx = (xmax-xmin)/n;
    xv = xmin:dx:xmax;
    x = round((x-xmin)/dx)+1;
    xmax = round((xmax-xmin)/dx)+1;
end   
z = zeros(length(xv)*size(x,2),1);
if nargin<3
    % make x one-dimensional by adding 1*xmax to the second column of x,
    % 2*xmax to the 3rd column of x,..., then sorting each column separately in
    % ascending order, then putting them all underneath 1st column -> save as num
    num = sort(x+xmax*ones(size(x,1),1)*(0:size(x,2)-1));
    num = num(:);
    % mark in 'z' which of the bins are occupied by at least one entry of x
    z(num) = 1;
    % diff([0;num;0])==0 is true if adjacent entries belong to the same bin
    % if several entries belong to the same bin, tmp==1 at the first, 
    % tmp==-1 at the last and tmp==0 everywhere else (tmp==1 and tmp==-1 have
    % same number of occurences!)
    tmp = diff(diff([0; num; 0])==0);
    ind = (1:length(num))'; % ascending numbers 1,2,3,... (auxiliary variable)
    % z now contains zeros and ones, bins that contain more than one entry
    % have to bet set to the right number: ind(tmp==-1)-ind(tmp==1)=no. of
    % repetitions of the bin-1 (but z was already 1 before, so it fits)
    z(num(tmp==1)) = z(num(tmp==1))-ind(tmp==1)+ind(tmp==-1);
else
    num = x+xmax*ones(size(x,1),1)*(0:size(x,2)-1);
    num = num(:);
    [num, ord] = sort(num);
    weight = weight(:);
    weight = weight(ord);
    tmp = diff([0; num])>0;
    z(num(tmp)) = weight(tmp);    
    tmp = diff(diff([0; num; 0])==0);
    ind = cumsum(weight);
    z(num(tmp==1)) = z(num(tmp==1)) + ind(tmp==-1)-ind(tmp==1); 
end    
z = reshape(z,length(xv),size(x,2));
if nargout==0
    bar(xv,z,1,'b'); 
    clear z xv
end
