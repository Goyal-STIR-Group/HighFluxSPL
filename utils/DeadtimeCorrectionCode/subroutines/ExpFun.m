function [err, c, zz, z] = ExpFun(p, t, y, pic, nrm, pos)

if length(t)<length(y)
    c = t;
    t = y(:);
    p = 1./p(:)';
    t = t(:);%-min(abs(t(:))); 
    zz = [ones(length(t),1) exp(-abs(t)*p)];
    if nargin>4 && ~isempty(nrm) 
        zz = zz./(ones(length(t),1)*sum(zz));
    end
    for j=1:size(c,2)
        err(:,j) = zz*c(:,j);
    end
else
    t = t(:);%-min(abs(t(:))); 
    p = 1./p(:)';
    [m, n] = size(y);
    if m<n y=y'; tmp=n; n=m; m=tmp; end
	t = t(isfinite(sum(y,2)));
    y = y(isfinite(sum(y,2)),:);
    [m, n] = size(y);
    zz = [ones(m,1) exp(-abs(t)*p)];
    if nargin>4 && ~isempty(nrm) 
        zz = zz./(ones(m,1)*sum(zz));
    end

    for j=1:n
        if nargin>5 && ~isempty(pos)
            c(:,j) = lsqnonneg(zz,y(:,j));
        else
            c(:,j) = zz\y(:,j);
        end
        z(:,j) = zz*c(:,j);
    end

    if nargin>3 && ~isempty(pic)
        if pic==1
            subplot(4,1,1:3)
            plot(t, y, 'ob', t, z, 'r'); 
            axis tight
            subplot(4,1,4)
            plot(t, (y-z)./sqrt(abs(z)),t,0*t)
            axis tight
            drawnow
        elseif pic==2
            semilogy(t, y, 'o', t, z); drawnow
        else
            semilogx(t, y, 'o', t, z); drawnow
        end
    end

    err = sum(sum((y-z).^2./abs(z)));
    %err = sum((y-z).^2);
    %ind = y>0;
    %err = sum(y(ind).*log(y(ind)./z(ind))-y(ind)+z(ind));
end

