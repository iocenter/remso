function [  shootingSol,Jacs,convergences ] = mrstSimulationStep( shootingVars,reservoirP,varargin)
%
%  MRST simulator function
%


opt = struct('shootingGuess',[],'force_step',false,'stop_if_not_converged',false);
opt = merge_options(opt, varargin{:});

[shootingSol.wellSols,shootingSol.ForwardStates,shootingSol.schedule,~,convergence,Jacs] = runScheduleADI(shootingVars.state0,...
                                                                                reservoirP.G,...
                                                                                reservoirP.rock,...
                                                                                reservoirP.system,...
                                                                                shootingVars.schedule,...
                                                                                'stop_if_not_converged', opt.stop_if_not_converged,...
                                                                                'initialGuess',opt.shootingGuess,...
                                                                                'force_step',opt.force_step);
                                                                            
                                                                            
convergences.residuals =  vertcat(convergence.residuals);
convergences.its =  sum(vertcat(convergence.its));
convergences.converged = all(vertcat(convergence.converged));
                                                                            
                                                                            
end

