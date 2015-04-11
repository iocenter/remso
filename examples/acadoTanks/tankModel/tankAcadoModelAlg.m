function [varargout] = tankAcadoModelAlg(xStart,u,dt,varargin)

varargout = cell(1,nargout);
%xp,dxpdu,dxpd0,convergence,shootingSol
opt = struct('gradients',false,'xLeftSeed',[],'vLeftSeed',[],'xRightSeeds',[],'guessV',[],'guessX',[],'uRightSeeds',[],'simVars',[],'p',0.4);
opt = merge_options(opt, varargin{:});

settings = ACADOintegratorsSettings;

settings.Model = 'tanks';
settings.Integrator = 'RK45';


settings.Tolerance = 1e-9;     % local error tolerance.
settings.AbsoluteTolerance = 1e-9;     % local error tolerance.
settings.MinimumStepSize = 1e-9;

settings.u = u;
settings.p = opt.p;

tStart = 0;
tEnd   = dt;

if opt.gradients
    
    if ~(size(opt.xLeftSeed,2)==0) && ~(size(opt.xRightSeeds,1)==0)
        error('not implemented')
    elseif (size(opt.xLeftSeed,2)==0) && (size(opt.xRightSeeds,1)==0)
%         settings.SensitivityMode = 'AD_FORWARD';  
%         settings.lambdaX         =  opt.xRightSeeds;
%         settings.lambdaU         =  opt.uRightSeeds;

        settings.SensitivityMode = 'AD_BACKWARD';
        settings.mu         =  eye( length(xStart) );
    elseif ~(size(opt.xLeftSeed,2)==0)
        settings.SensitivityMode = 'AD_BACKWARD';
        settings.mu         =  eye( length(xStart) );
    elseif ~(size(opt.xRightSeeds,1)==0)
%         settings.SensitivityMode = 'AD_FORWARD';  
%         settings.lambdaX         =  opt.xRightSeeds;
%         settings.lambdaU         =  opt.uRightSeeds;

        % possible bug in acado multiply afterwards!
        settings.SensitivityMode = 'AD_BACKWARD';
        settings.mu         =  eye( length(xStart) );

    else
        error('what')
    end
end


[ xEnd, outputsB ] = ACADOintegrators( settings,xStart,tStart,tEnd );

varargout{1} = xEnd;
varargout{2} = xEnd(1)-xStart(1);
varargout{3} = [];




if opt.gradients

    vJx = outputsB.Jx(1,:);
    vJx(1,1) = vJx(1,1) - 1; 
    vJu = outputsB.Ju(1,:);
    
    if ~(size(opt.xLeftSeed,2)==0) && ~(size(opt.xRightSeeds,1)==0)
        error('not implemented')
    elseif (size(opt.xLeftSeed,2)==0) && (size(opt.xRightSeeds,1)==0)
        % should be! settings.SensitivityMode = 'AD_FORWARD';  
        
        Jac.xJx = outputsB.Jx;
        Jac.xJu = outputsB.Ju;
        Jac.vJx = vJx;
        Jac.vJu = vJu;
        
    elseif ~(size(opt.xLeftSeed,2)==0)
        Jac.Jx = opt.xLeftSeed*outputsB.Jx + opt.vLeftSeed * vJx;
        Jac.Ju = opt.xLeftSeed*outputsB.Ju + opt.vLeftSeed * vJu;
        
    elseif ~(size(opt.xRightSeeds,1)==0)
                % should be! settings.SensitivityMode = 'AD_FORWARD';  

        Jac.xJ = outputsB.Jx*opt.xRightSeeds + outputsB.Ju*opt.uRightSeeds;
        Jac.vJ = vJx*opt.xRightSeeds + vJu*opt.uRightSeeds;
    else
        error('what')
    end
    varargout{3} = Jac;
end
if outputsB.Status == 0
    varargout{4} = struct('converged',true);
else
    varargout{4} = struct('converged',false);    
end
 varargout{5} = [];    


