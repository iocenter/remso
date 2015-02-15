function [yi, dyidxi, dyidvi] = interp2DPVT(T, xi, vi, reginx)
compDer = (nargout>1);
nreg = numel(reginx);

%{
Changes by codas:

This are changes suggested by stein to mrst2014a but that did not appear in mrst2014b
why?
%}


w   = zeros(size(xi));
a   = zeros(size(xi));
dx  = zeros(size(xi));
yil = zeros(size(xi));
yir = zeros(size(xi));
if compDer
    dyidxil = zeros(size(xi));
    dyidxir = zeros(size(xi));
    dwdvi   = zeros(size(xi));
end

for k = 1:nreg
    v = T{k}.key;  pos = T{k}.pos;
    x = T{k}.data(pos(1:end-1),1);
    lims = v; lims(1)=-inf; lims(end)=inf;
    vi_k = vi(reginx{k});
    [bin,bin] = histc(vi_k, lims);                                     %#ok
    w(reginx{k}) = (vi_k-v(bin))./(v(bin+1)-v(bin));
    dx(reginx{k}) = (x(bin+1)-x(bin));
    a(reginx{k}) = (x(bin+1)-x(bin))./(v(bin+1)-v(bin));
    if compDer
        dwdvi(reginx{k}) = 1./(v(bin+1)-v(bin));
    end
    for tn = 1:numel(v)
        lns = pos(tn):(pos(tn+1)-1);
        tab = T{k}.data(lns,:);

        ixl = (bin==(tn-1));
        ixr = (bin==tn);
        if ~(ischar(reginx{k}) && reginx{k} == ':')
%             ixl = ixl.*reginx{k};
%             ixr = ixr.*reginx{k};
            ixl = reginx{k}(ixl);
            ixr = reginx{k}(ixr);        
        end
        % shift values such that intp is linear along (x,v)
        xils = xi(ixl) + (1-w(ixl)).*dx(ixl);
        xirs = xi(ixr) - w(ixr).*dx(ixr);
        
        yil(ixl) = interpTable(tab(:,1), tab(:,2), xils);
        yir(ixr) = interpTable(tab(:,1), tab(:,2), xirs);
        if compDer
            dyidxil(ixl) = dinterpTable(tab(:,1), tab(:,2), xils);
            dyidxir(ixr) = dinterpTable(tab(:,1), tab(:,2), xirs);
        end
    end
end
yi = yil.*w + yir.*(1-w);
if compDer
    dyidxi = dyidxil.*w + dyidxir.*(1-w);
    dyidvi = (yil-yir).*dwdvi - dyidxi.*a;
end
end
