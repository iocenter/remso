function [ f,df,vars,simVars,debugInfo ] = lineFunctionWrapper(stepL,x0,v0,u0,s0,dx,dv,du,ds,simF,target,meritF,jobSchedule,varargin)
%
% Calculate the value of the merit fuction on the line
%
% SYNOPSIS:
%  [f,df,vars,simVars,debugInfo] = lineFunctionWrapper(stepL,x0,v0,u,dx,dv,du,simF,target,meritF,jobSchedule, obj)
%  [f,df,vars,simVars,debugInfo] = lineFunctionWrapper(stepL,x0,v0,u,dx,dv,du,simF,target,meritF,jobSchedule, 'pn', pv, ...)
% PARAMETERS:
%
%  stepL - value from 0 to 1, point on the line to evaluate
%
%  (x0,v0,u0) - point at the beguining of the line
%
%  (dx,dv,du) - line vector
%
%  simF - multiple shooting simulator function
%
%  target - target function representing the objective
%
%  meritF - merit function
%
%  jobSchedule - Object with parallel workers information
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%  gradients - if true compute the gradient on the line
%
%  fwd - if true select forward mode AD.
%
%  simVars - simVars for the current point.
%
%  (xd0,vd0) - (xs-x,vs-x) at the beguining of the line.
%
%  (xs0,vs0) - (xs,vs) at the beguining of the line. 
%
%  plotFunc, plot, debug, set this options for debugging
%
% RETURNS:   
%
%   f - line function value
%
%   fd - line function gradients.
% 
%   vars - (x,v,u) at the end of the lien
%
%   simVars - simulator variables at the evaluated point
%
% SEE ALSO:
%
%

opt = struct('gradients',true,'fwd',true,'simVars',[],'xd0',[],'vd0',[],'sd0',[],'xs0',[],'vs0',[],'s20',[],'plotFunc',[],'xi',0,'plot',false,'debug',false);
opt = merge_options(opt, varargin{:});

imMaster = jobSchedule.imMaster;

simVars = opt.simVars;
xi = opt.xi;
if opt.fwd && opt.gradients
    
    %  At the step 0 there is no need of simulation, we have all required
    %  information
    if stepL == 0
        
        x = x0;
        v = v0;
        s = s0;
        u = u0;
        
        xs = opt.xs0;
        vs = opt.vs0;
        s2 = opt.s20;
        
        xd0 = opt.xd0;
        vd0 = opt.vd0;
        sd0 = opt.sd0;
        
        %spmd
        xJ = diffJac(dx,xd0,xi);
        vJ = diffJac(dv,vd0,xi);
        %end
        if imMaster
        sJ = diffJacS(ds,sd0,xi);
        end
        
        usliced = [];
        
    else
        % simulate the new point
        
        % generate the point and a simulation guess
        xs0 = opt.xs0;
        vs0 = opt.vs0;
        s20 = opt.s20;

        u = cellfun(@(z,dz)z+stepL*dz,u0,du,'UniformOutput',false);
        
            %spmd
            x = walkLine(x0,dx,stepL);
            v = walkLine(v0,dv,stepL);
            %end
            
            s = walkLineS(s0,ds,stepL);

            if isempty(xs0)
                guessX = x;
            else
                %spmd
                guessX = buildGuess(xs0,x0,dx,stepL);
                %end
            end
            if isempty(vs0)
                guessV = v;
            else
                %spmd
                guessV = buildGuess(vs0,v0,dv,stepL);
                %end
            end
            
        
      
        
        % Multiple shooting simulation call!
        [xs,vs,s2,JacRes,convergence,simVars,usliced] = simF(x,u,v,'gradients',true,'guessX',guessX,'guessV',guessV,'xRightSeed',dx,'vRightSeed',dv,'uRightSeed',du,'simVars',simVars);
        
        xJ = JacRes.xJ;
        vJ = JacRes.vJ;
        sJ = JacRes.sJ;
          
    end
    

    
        %spmd
        eX = MINUS(xs,x);
        jeX = MINUS(xJ,dx);
        
        eV = MINUS(vs,v);
        jeV = MINUS(vJ,dv);
        %end
        if imMaster
        eS = minus(s2,s);
        jeS = minus(sJ,ds);
        else
        eS = nan;
        jeS =nan;
        end
        
    
    % objective function evalutation
    if imMaster
    [tarL,JacTar] = target(s,u,'gradients',true,'sRightSeed',ds,'uRightSeed',du);
    else
	tarL = nan;
    JacTar.J = nan; % force call of gradients
    end
    % merit function evalutaion
	[f,JacF,debugInfo] = meritF(tarL,{eX;eV},eS,'gradients',true,'fRightSeeds',JacTar.J,'dERightSeeds',{jeX;jeV},'dSRightSeeds',jeS,'debug',opt.debug);

    
    df = JacF.J;
    
	fdf = NMPI_Bcast([f,df],2,jobSchedule.Master_rank,jobSchedule.my_rank);
    f =fdf(1);
    df =fdf(2);
    
elseif opt.gradients
    error('not Implemented')
else
    error('not Implemented')
end

vars.x=x;
vars.u=u;
vars.v=v;
vars.s=s;


vars.xs = xs;
vars.vs = vs;
vars.s2 = s2;
vars.usliced = usliced;


if ~isempty(opt.plotFunc) && opt.plot
        opt.plotFunc(vars.x,...
            vars.u,...
            vars.v,...
            vars.s,...
            eX);
end

if isempty(simVars);
    itsMean = nan;
else
    [its] = getSolverIterations(simVars);
	sumIts = sum(its);
    n = numel(its);
	sumItsN = gopMPI('+',[sumIts,n],jobSchedule);
    itsMean = sumItsN(1)/sumItsN(2);
end
debugInfo.itsMean = itsMean;

end


function [x] = walkLine(x0,dx,stepL)
    f = @(z,dz)walkLineS(z,dz,stepL);   
    x = cellfun(@(zr,dzr)cellfun(f,zr,dzr,'UniformOutput',false),x0,dx,'UniformOutput',false);
end
function [x] = walkLineS(z,dz,stepL)
    x = z+stepL*dz;
end


function [guessX] = buildGuess(xs0,x0,dx,stepL)
    f = @(zs0,z0,dz)buildGuessS(zs0,z0,dz,stepL);
	guessX = cellfun(@(zs0r,z0r,dzr)cellfun(f,zs0r,z0r,dzr,'UniformOutput',false),xs0,x0,dx,'UniformOutput',false);
end
function [guessX] = buildGuessS(zs0,z0,dz,stepL)
	guessX = zs0*(1-stepL)+(z0+dz)*stepL;
end


function [zJ] = diffJac(dz,zd0,xi)
    f = @(dz,zd0)diffJacS(dz,zd0,xi);
	zJ = cellfun(@(x1r,x2r)cellfun(f,x1r,x2r,'UniformOutput',false),dz,zd0,'UniformOutput',false);
end
function [zJ] = diffJacS(dz,zd0,xi)
	zJ = (dz-(1-xi)*zd0);
end


function [zJ] = MINUS(z1,z2)
    f = @minus;
	zJ = cellfun(@(z1r,z2r)cellfun(f,z1r,z2r,'UniformOutput',false),z1,z2,'UniformOutput',false);
end

