function [  shootingSol,Jac,convergence ] = mrstSimulationStep( shootingVars,reservoirP,varargin)
%
%  MRST simulator function
%


opt = struct('shootingGuess',[],'stop_if_not_converged', false);
opt = merge_options(opt, varargin{:});

[shootingSol.wellSols,shootingSol.ForwardStates,shootingSol.schedule,~,convergence,Jac] = runScheduleADI(shootingVars.state0,...
                                                                                reservoirP.G,...
                                                                                reservoirP.rock,...
                                                                                reservoirP.system,...
                                                                                shootingVars.schedule,...
                                                                                'stop_if_not_converged', opt.stop_if_not_converged, ...
                                                                                'initialGuess',opt.shootingGuess);


end

