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

opt = struct('m',6,'scale',false,'dF',0.2,'epsd',1e-5,'condT',1e9,'debug',true,'it',0);
opt = merge_options(opt, varargin{:});

if opt.debug
	fid = fopen('logBFGS.txt','a');   
end

stepNorm = norm(du);
skipping = stepNorm < opt.epsd;
if skipping;
    if opt.debug
        fprintf(fid,'%3.d BFGS not updated: too short step %e \n',opt.it,stepNorm);
    end
    return;
end


if opt.m <= 0  %%% do not save any curvature!
    
elseif size(S,2) >= opt.m
    % pick the last opt.m-1 elements and add the currect curvature at the end
    mL = size(S,2);
    mI = mL - opt.m+2;
    S = [S(:,mI:mL),du];
    Y = [Y(:,mI:mL),yG'];
else
    S = [S,du];
    Y = [Y,yG'];
end

if isempty(M)
    [M,condM,scaledOk] = initializeBFGS(nru,opt.scale,yG,du,opt.condT);
    if opt.debug && ~scaledOk
        fprintf(fid,'%3.d BFGS scaling failed during restart\n',opt.it);
    end
end

[ M,skipping,damping,minEig ] = dampedBfgsUpdate(M,yG,du,'dF',opt.dF,'epsd',opt.epsd);

if opt.debug && skipping
    % if we skipped it is because of negative curvature, which is imposible!
    fprintf(fid,'%3.d Check the BFGS code min(eig(Q*))= %1.1e',opt.it,minEig);
end
if opt.debug && damping
	fprintf(fid,'%3.d Damped BFGS update\n',opt.it);
end

%% Restart the Hessian approximation!
% skipping == true if there is a negative eigenvalue!
condM = cond(M);
if condM > opt.condT
    if opt.debug && condM > opt.condT
        fprintf(fid,'%3.d BFGS Restart, |Q*| %1.1e, ',opt.it,condM);
    end
    
    [M,condM,scaledOk] = initializeBFGS(nru,opt.scale,yG,du,opt.condT);
    if opt.debug && ~scaledOk
        fprintf(fid,'%3.d BFGS scaling failed during restart\n',opt.it);
    end
    
    skippingCounter = 0;
    dampingCounter = 0;
    for k = size(S,2):-1:1
        [ MT,skipping,damping,minEig ] = dampedBfgsUpdate(M,(Y(:,k))',S(:,k),'dF',opt.dF,'epsd',opt.epsd);
        
        % Apply the update if the resulting condition do not exceed the
        % threshold
        if skipping
            skippingCounter = skippingCounter + 1;
        end
        if damping
            dampingCounter = dampingCounter + 1;
        end
        
        condMT = cond(MT);
        if condMT < opt.condT
            condM = condMT;
            M = MT;
        end
    end
    
    if opt.debug
        fprintf(fid,'|Q| %1.1e dampCount %d skipCount %d\n',condM,dampingCounter,skippingCounter);
    end
    
end

if opt.debug
    fclose(fid);
end

end




function [M,condM,scaledOk] = initializeBFGS(nru,scale,yG,du,condT)
scaledOk = true;
if scale
    sTy = dot(yG,du);
    if sTy > 0
        M = eye(nru)*(yG'*yG)/sTy;
        condM = cond(M);
        if condM > condT
            scaledOk = false;
            M = eye(nru);
            condM = 1;
        end
    else
        scaledOk = false;
        M = eye(nru);
        condM = 1;
    end
else
    M = eye(nru);
    condM = 1;
end
end