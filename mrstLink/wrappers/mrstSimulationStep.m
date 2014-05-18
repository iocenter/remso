function [  shootingSol,Jac,convergence ] = mrstSimulationStep( shootingVars,reservoirP,varargin)
%
%  MRST simulator function
%


opt = struct('shootingGuess',[]);
opt = merge_options(opt, varargin{:});

[shootingSol.wellSols,shootingSol.ForwardStates,~,convergence,Jac] = runScheduleADI(shootingVars.state0,...
                                                                                reservoirP.G,...
                                                                                reservoirP.rock,...
                                                                                reservoirP.system,...
                                                                                shootingVars.schedule,...
                                                                                'initialGuess',opt.shootingGuess);


end

