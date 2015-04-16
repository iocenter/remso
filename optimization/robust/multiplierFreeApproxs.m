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


gbarR = gbar.Js*as + innerProd(gbar.Jx,ax) + innerProd(gbar.Jv,av);



d = [vertcat(xd{:});vertcat(vd{:});{sd}];
errorSum = sum(cellfun(@(di)sum(abs(di)),d));
crossProduct = sum(cellfun(@(wi,dui)wi*dui,w,du'))/(1-xi);






end


function ip = innerProd(J,a)

ip = sum(cellfun(@innerProdS,J,a'));

end

function ip = innerProdS(J,a)

ip = cell2mat(J)*cell2mat(a);

end