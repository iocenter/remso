function [ error ] = testRunAdjointADI( initState, G, rock, fluid, schedule, system,xScale,vScale,cellControlScales,varargin)
%  
%  TEST modified runAdjointADI against finite differences. Try it in a small
%  model, it can be expensive.  
%
%  Recommended to run when updating to a new MRST release.  Observe that
%  the gradients given by runAdjointADI are TRANSPOSE with respect to the
%  gradients in REMSO.

opt = struct('pert', 1e-4, 'finalVars', true);
opt = merge_options(opt, varargin{:});

nw = numel(schedule.control(1).W);
nvw = nw*3;

if opt.finalVars
    finalTime = schedule.time+sum(schedule.step.val);
    f = @(k,wellSols,states,scheduleOut,gradFlag) finalStepVars(k,states{end},wellSols{end},scheduleOut,finalTime,'ComputePartials',gradFlag,'xvScale',[xScale;vScale(1:nvw)]);
else
    f = @(k,wellSols,states,scheduleOut,gradFlag) unPack(NPVOW(G, wellSols, scheduleOut,  'ComputePartials',gradFlag,'tStep',k));
end

pert = opt.pert;

x0 = stateMrst2stateVector( initState,'xScale',xScale );

uScale = cellControls2Controls(cellControlScales);
u  = schedule2Controls( schedule,'uScale',uScale);


[obj,jx,ju] = fRemso(x0,u,xScale,uScale,f, G, rock, fluid, schedule, system,true);


fx = @(x0)fRemso(x0,u,xScale,uScale,f, G, rock, fluid, schedule, system,false);
fu = @(u) fRemso(x0,u,xScale,uScale,f, G, rock, fluid, schedule, system,false);

[dfdu] = calcPertGrad(fu,u,pert);

[dfdx] = calcPertGrad(fx,x0,pert);


error = max(max(max(abs(dfdx-jx))),max(max(abs(ju-dfdu))));

end


function [obj,jx,ju] = fRemso(x0,u,xScale,uScale,f, G, rock, fluid, schedule, system,gradFlag)


[ initState ] = stateVector2stateMrst( x0,'xScale',xScale);
[ schedule ] = controls2Schedule( u,schedule,'uScale',uScale);

[obj,grad] = fMRST(f, initState, G, rock, fluid, schedule, system,gradFlag);

jx = [];
ju = [];

if gradFlag
    jx = bsxfun(@times,grad{1}',xScale');
    ju = cell2mat(grad(2:end)');
    ju = bsxfun(@times,ju',uScale');
end

end



function [obj,grad] = fMRST(f, initState, G, rock, fluid, schedule, system,gradFlag)

[wellSols,states,scheduleOut,iter,convergence]= runScheduleADI(initState, G, rock, system, schedule);

objective = @(k) f(k,wellSols,states,scheduleOut,gradFlag);


obj = objective(1);
for k = 2:numel(schedule.step.val)
    obj = obj + objective(k) ;
end

grad = [];
if gradFlag
    grad = runAdjointADI(G, rock, fluid, scheduleOut, objective, system,'ForwardStates',[{initState};states],....
        'initialConditionSens', true);
end

end

function [x] = unPack(objk)

if iscell(objk)
    x = objk{:};
end

end


function [dfdx] = calcPertGrad(f,xb,pert)

nx = numel(xb);
fb = f(xb);

dfdx = zeros(numel(fb),nx);

parfor j = 1:nx
    xp = xb;
    xp(j) = xp(j) + pert;
    
    fp = f(xp);
    
    dfdx(:,j) = (fp-fb)/pert;
    
end

end