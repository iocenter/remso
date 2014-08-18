function [  shootingSol,Jacs,convergence ] = mrstSimulationStep( shootingVars,reservoirP,varargin)
%
%  MRST simulator function
%


opt = struct('shootingGuess',[],'force_step',false);
opt = merge_options(opt, varargin{:});

[shootingSol.wellSols,shootingSol.ForwardStates,shootingSol.schedule,~,convergence,Jacs] = runScheduleADI(shootingVars.state0,...
                                                                                reservoirP.G,...
                                                                                reservoirP.rock,...
                                                                                reservoirP.system,...
                                                                                shootingVars.schedule,...
                                                                                'initialGuess',opt.shootingGuess,'force_step',opt.force_step);


end

