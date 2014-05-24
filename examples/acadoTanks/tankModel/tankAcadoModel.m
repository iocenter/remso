function [varargout] = tankAcadoModel(xStart,u,dt,varargin)

varargout = cell(1,nargout);
%xp,dxpdu,dxpd0,convergence,shootingSol
opt = struct('gradients',false,'xLeftSeed',[],'vLeftSeed',[],'xRightSeeds',[],'guessV',[],'guessX',[],'uRightSeeds',[],'simVars',[]);
opt = merge_options(opt, varargin{:});

settings = ACADOintegratorsSettings;

settings.Model = 'tanks';
settings.Integrator = 'RK45';


settings.Tolerance = 1e-11;     % local error tolerance.
settings.AbsoluteTolerance = 1e-11;     % local error tolerance.
settings.MinimumStepSize = 1e-9;

settings.u = u;


tStart = 0;
tEnd   = dt;


if opt.gradients
    
    if ~isempty(opt.xLeftSeed) && ~isempty(opt.xRightSeeds)
        error('not implemented')
    elseif isempty(opt.xLeftSeed) && isempty(opt.xRightSeeds)
%         settings.SensitivityMode = 'AD_FORWARD';  
%         settings.lambdaX         =  opt.xRightSeeds;
%         settings.lambdaU         =  opt.uRightSeeds;

        settings.SensitivityMode = 'AD_BACKWARD';
        settings.mu         =  eye( length(xStart) );
    elseif ~isempty(opt.xLeftSeed)
        settings.SensitivityMode = 'AD_BACKWARD';
        settings.mu         =  opt.xLeftSeed;
    elseif ~isempty(opt.xRightSeeds)
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
    
    %settings.SensitivityMode = 'AD_FORWARD';
    %settings.lambdaU         =  ones(length(xStart) );  % forward seed
    
    
    [ xEnd, outputsB ] = ACADOintegrators( settings,xStart,tStart,tEnd );
    
    varargout{1} = xEnd;
    varargout{2} = [];
    varargout{3} = [];
    
if opt.gradients
    
    if ~isempty(opt.xLeftSeed) && ~isempty(opt.xRightSeeds)
        error('not implemented')
    elseif isempty(opt.xLeftSeed) && isempty(opt.xRightSeeds)
        % should be! settings.SensitivityMode = 'AD_FORWARD';  
        
        Jac.xJx = outputsB.Jx;
        Jac.xJu = outputsB.Ju;
    elseif ~isempty(opt.xLeftSeed)
        Jac.Jx = outputsB.Jx;
        Jac.Ju = outputsB.Ju;
        
    elseif ~isempty(opt.xRightSeeds)
                % should be! settings.SensitivityMode = 'AD_FORWARD';  

        Jac.xJ = outputsB.Jx*opt.xRightSeeds + outputsB.Ju*opt.uRightSeeds;
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

    