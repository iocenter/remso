function [gbarR,errorSum,crossProduct] = multiplierFreeApproxs(gbar,ax,av,as,xd,vd,sd,w,du,xi)
%
%
% This values will generate a mu that induces a descent direction for zeta
% 
% gbarR = gbar Y p_y
% errorSum = |c|_1
% crossProduct = w'*p_z/(1-\xi)
%


% In robust optimization these are zero

%gbarR =  sum(cellfun(@(gbari,dzi)gbari*dzi,gbar.Jx',ax));
%au = 0 !
%if withAlgs    
%    gbarR = gbarR + sum(cellfun(@(gbari,dzi)gbari*dzi,gbar.Jv',av));
%end

gbarJx = gbar.Jx;
gbarJv = gbar.Jv;
spmd
gbarR = innerProd(gbarJx,ax) + innerProd(gbarJv,av);
gbarR = gop(@plus,gbarR);
end
gbarR = gbar.Js*as + gbarR{1};

spmd
    errorSum = errorSumZ(xd) + errorSumZ(vd);
    errorSum = gop(@plus,errorSum);
end
errorSum = errorSum{1};

errorSum = errorSum + sum(abs(sd));
crossProduct = sum(cellfun(@(wi,dui)wi*dui,w,du'))/(1-xi);






end


function ip = innerProd(J,a)

ip = sum(cellfun(@innerProdS,J,a));

end

function ip = innerProdS(J,a)

ip = cell2mat(J)*cell2mat(a);

end

function eS = errorSumZ(z)
if ~isempty(z)
    eS = sum([cellfun(@(di)sum(abs(di)),vertcat(z{:}));0]);
else
    eS = 0;
end
end