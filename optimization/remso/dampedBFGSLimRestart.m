function [  M,S,Y, skipping ] = dampedBFGSLimRestart(M,yG,du,nru,S,Y,varargin )
% Dampeg BFGS Hessian approximation with Limited memory restart
%
% SYNOPSIS:
%  [M,S,Y, skipping] = dampedBFGSLimRestart(M,yG,du,nru,S,Y)
%  [M,S,Y, skipping] = dampedBFGSLimRestart(M,yG,du,nru,S,Y, 'pn', pv, ...)
% PARAMETERS:
%
%   M - Current hessian approximation
%
%   yG - Lagrangian gradient difference.
%
%   du - variables step difference
%
%   nru - number of variables
%
%   S - 'du' memory
%
%   Y - 'yG' memory
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%   m - Number of previous steps to keep in memory for restart.
%
%   scale - Reinitialize and scale the hessian (called in the first step)
%
%   dF - Damping factor according to M.J.D Powell
%
%   epsd - variables minimum step threshold.
%
%   condT - Hessian approximation condition threshold.
%
%   debug - if true, print debug information
%
%   it - it number for debug
%
% RETURNS:
%
%   M - Updated Hessian approximation
%   
%   S - Updated 'du' memory 
%
%   Y - Updated 'Y' memory 
%
% SEE ALSO:
%
%

opt = struct('m',6,'scale',false,'dF',0.2,'epsd',1e-5,'condT',1e12,'debug',true,'it',0);
opt = merge_options(opt, varargin{:});

if opt.debug
    if opt.it <= 1
        fid = fopen('logBFGS.txt','w');
    else
        fid = fopen('logBFGS.txt','a');
        
    end
else
    
end

stepNorm = norm(du);
skipping = stepNorm < opt.epsd;
if skipping;
    if opt.debug
        fprintf(fid,'%d BFGS not updated: too short step %e \n',opt.it,stepNorm);
    end
    return;
end

if size(S,2) >= opt.m
    S = [S(:,2:end),du];
    Y = [Y(:,2:end),yG'];
else
    S = [S,du];
    Y = [Y,yG'];
end


[ M,skipping,minEig ] = dampedBfgsUpdate(M,yG,du,nru,'scale',opt.scale,'dF',opt.dF,'epsd',opt.epsd);


%% Restart the Hessian approximation!
% skipping == true if there is a negative eigenvalue!
condM = cond(M);
if skipping  || condM > opt.condT
    if opt.debug && condM > opt.condT
        fprintf(fid,'%2.d BFGS Restart, |Q*| %1.1e, ',opt.it,condM);
    end
    if opt.debug && skipping
        fprintf(fid,'%2.d BFGS Restart, min(eig(Q*))= %1.1e',opt.it,minEig);
    end
    
    M = eye(nru);
    condM = 1;
    scale = true;
    for k = size(S,2):-1:1
        [MT] = dampedBfgsUpdate(M,(Y(:,k))',S(:,k),nru,'scale',scale,'dF',opt.dF,'epsd',opt.epsd);
        
        % Apply the update if the resulting condition do not exceed the
        % threshold
        condMT = cond(MT);
        if condMT < opt.condT
            condM = condMT;
            M = MT;
        end
        scale = false;
    end
    
    if opt.debug
        fprintf(fid,'|Q| %1.1e\n',condM);
    end
    
end

if opt.debug
    fclose(fid);
end

end

