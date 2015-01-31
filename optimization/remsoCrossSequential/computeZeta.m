function [ zeta ] = computeZeta( gTZ,B,w )

% L. Biegler, C. Schmid, and D. Ternet, 'A multiplier-free, reduced Hessian method for process optimization' in Large-Scale
% Optimization with Applications, ser. The IMA Volumes in Mathematics and its Applications, L. T. Biegler, T. Coleman,
% A. R. Conn, and F. Santosa, Eds. Springer New York, 1997, vol. 93, pp. 101–127.


if iscell(gTZ)
    gTZ = cell2mat(gTZ);
end
if iscell(w)
    w = cell2mat(w);

gTZBinvT = B\gTZ';

gTZBinvw =  w*gTZBinvT;

if gTZBinvw >= 0
    zeta = 1;
else
    zeta = min(-0.1*gTZ*gTZBinvT/gTZBinvw,1);
end

end

