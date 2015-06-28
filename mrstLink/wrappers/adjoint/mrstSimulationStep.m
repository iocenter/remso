function [  shootingSol,Jacs,convergences ] = mrstSimulationStep( shootingVars,reservoirP,varargin)
%
%  MRST simulator function
%

opt = struct('shootingGuess',[]);
opt = merge_options(opt, varargin{:});

[simRes,reports] = runSchedule(shootingVars.state0,...
    reservoirP.G,...
    reservoirP.S,...
    reservoirP.W,...
    reservoirP.rock,...
    reservoirP.fluid,...
    shootingVars.schedule,...
    'init_state',opt.shootingGuess,...
    'gravityOff',true,...
    'VerboseLevel', 0);


shootingSol.ForwardStates = simRes(2:end);  % remove initial condition
shootingSol.schedule = shootingVars.schedule;


convergences.residuals =  vertcat(reports.residual);
convergences.its =  vertcat(reports.iterations);
convergences.converged = all(vertcat(reports.success));

Jacs = [];

end

