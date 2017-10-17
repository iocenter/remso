function [  shootingSol,Jacs,convergences ] = mrstSimulationStep( shootingVars,reservoirP,varargin)
%
%  MRST simulator function
%
%   TODO:
%  What about the Linear and nonlinear solvers??? where to place them?

opt = struct('shootingGuess',[]);
opt = merge_options(opt, varargin{:});


[shootingSol.wellSols, shootingSol.ForwardStates, schedulereport] = ...
    simulateScheduleAD(shootingVars.state0,...
                       reservoirP.model,...
                       shootingVars.schedule,...
                       'OutputMinisteps', true,...
                       'initialGuess',opt.shootingGuess);
                   
[shootingSol.schedule] = convertReportToSchedule(schedulereport, shootingVars.schedule);
                   
Jacs = [];
                                                                            
convergences.its =  sum(schedulereport.Iterations);
convergences.converged = all(schedulereport.Converged);
                                                                            
end

