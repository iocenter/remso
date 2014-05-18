function [ f,df,vars,simVars,debugInfo ] = lineFunctionWrapper(stepL,x0,v0,u0,dx,dv,du,simF,target,meritF,varargin)
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

opt = struct('gradients',true,'fwd',true,'simVars',[],'xd0',[],'vd0',[],'xs0',[],'vs0',[],'plotFunc',[],'plot',false,'debug',false);
opt = merge_options(opt, varargin{:});


withAlgs = ~isempty(v0{1});
simVars = opt.simVars;
if opt.fwd && opt.gradients
    
    %  At the step 0 there is no need of simulation, we have all required
    %  information
    if stepL == 0
        
        x = x0;
        if withAlgs
            v = v0;
        else
            v = [];
        end
        
        u = u0;
        
        xs = opt.xs0;
        
        
        xd0 = opt.xd0;
        vd0 = opt.vd0;
        vs = opt.vs0;
        
        
            xJ = cellfun(@minus,dx,xd0,'UniformOutput',false);
            
            if withAlgs
                vJ = cellfun(@minus,dv,vd0,'UniformOutput',false);
            else
                vs =[];
                vJ = [];
            end
            
            
        
        
        usliced = [];
        
    else
        % simulate the new point
        
        % generate the point and a simulation guess
        xs0 = opt.xs0;
        vs0 = opt.vs0;
        u = cellfun(@(z,dz)z+stepL*dz,u0,du,'UniformOutput',false);
        
            x = walkLine(x0,dx,stepL);
            if withAlgs
                v = walkLine(v0,dv,stepL);
            else
                v = [];
            end
            if isempty(xs0)
                guessX = x;
            else
                [guessX] = buildGuess(xs0,x0,dx,stepL);
            end
            if isempty(vs0)
                guessV = v;
            elseif withAlgs
                [guessV] = buildGuess(vs0,v0,dv,stepL);
            else
                guessV = [];
            end
            

        
      
        
        % Multiple shooting simulation call!
        [xs,vs,JacRes,convergence,simVars,usliced] = simF(x,u,'gradients',true,'guessX',guessX,'guessV',guessV,'xRightSeed',dx,'uRightSeed',du,'simVars',simVars);
        
        xJ = JacRes.xJ;
        vJ = JacRes.vJ;
          
    end
    

    
    
        eX = cellfun(@minus,xs,x,'UniformOutput',false);
        jeX = cellfun(@minus,xJ,dx,'UniformOutput',false);
        
        if withAlgs
            eV = cellfun(@minus,vs,v,'UniformOutput',false);
            jeV = cellfun(@minus,vJ,dv,'UniformOutput',false);
        else
            eV = [];
            jeV = [];
        end
        if ~withAlgs
            vs =[];
            vJ = [];
        end   
    
    % objective function evalutation
    [tarL,JacTar] = target(xs,u,vs,'gradients',true,'xRightSeeds',xJ,'uRightSeeds',du,'vRightSeeds',vJ);
    
    % merit function evalutaion
    [f,JacF,debugInfo] = meritF(tarL,{eX;eV},{xs;vs},'gradients',true,'fRightSeeds',JacTar.J,'dERightSeeds',{jeX;jeV},'bERightSeeds',{xJ;vJ},'debug',opt.debug);
    
    df = JacF.J;
    

    
elseif opt.gradients
    error('not Implemented')
else
    error('not Implemented')
end

vars.x=x;
vars.u=u;
vars.v=v;

vars.xs = xs;
vars.vs = vs;
vars.usliced = usliced;


if ~isempty(opt.plotFunc) && opt.plot
        opt.plotFunc(vars.x,...
            vars.u,...
            vars.v,...
            eX);
end

if isempty(simVars);
    itsMean = nan;
else
    [itsMean] = mean(getSolverIterations(simVars));
end
debugInfo.itsMean = itsMean;

end


function [x] = walkLine(x0,dx,stepL)
    x = cellfun(@(z,dz)z+stepL*dz,x0,dx,'UniformOutput',false);
end

function [guessX] = buildGuess(xs0,x0,dx,stepL)
	guessX = cellfun(@(zs0,z0,dz)zs0*(1-stepL)+(z0+dz)*stepL,xs0,x0,dx,'UniformOutput',false);
end

