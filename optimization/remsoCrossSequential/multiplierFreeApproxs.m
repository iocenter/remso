function [gbarR,errorSum,crossProduct] = multiplierFreeApproxs(gbar,ax,av,xd,vd,w,du,xi,withAlgs)
%
%
% This values will generate a mu that induces a descent direction for zeta
% 
% gbarR = gbar Y p_y
% errorSum = |c|_1
% crossProduct = w'*p_z/(1-\xi)
%




gbarR =  sum(cellfun(@(gbari,dzi)gbari*dzi,gbar.x',ax));
%au = 0 !
if withAlgs    
    gbarR = gbarR + sum(cellfun(@(gbari,dzi)gbari*dzi,gbar.v',av));
end




d = [xd;vd];
errorSum = sum(cellfun(@(di)sum(abs(di)),d));
crossProduct = sum(cellfun(@(wi,dui)wi.*dui,w,du'))/(1-xi);






end
