function [s,Jac,o ] = realization2s(x,u,v,sss,varargin )
% TODO: rename this function
% compute risk measure s given the variables v of each realization

opt1 = struct('partials',false,'vRightSeed',[],'xRightSeed',[],'uRightSeed',[]);
[opt1,others] = merge_options(opt1, varargin{:});

opt2 = struct('partials',opt1.partials,'leftSeed',[]);  % bugfix for the other codes
opt2 = merge_options(opt2, others{:});

imMaster = sss.jobSchedule.imMaster;

s = nan;
Jac = [];
if ~opt1.partials
    
	[o] = realizationOutput(x,u,v,sss,'partials',false);
    if imMaster
	[s] = outputRisks(o,'eta',sss.eta,'partials',false);    
    end
else
    if size(opt2.leftSeed,2)~=0  %% backward mode
        
        [o] = realizationOutput(x,u,v,sss,'partials',false);
        if imMaster
        [s,sJo] = outputRisks(o,'eta',sss.eta,'partials',true,'leftSeed',opt2.leftSeed);   
        else
        sJo.Jo = cell(0);
        end
        [o,Jac] = realizationOutput(x,u,v,sss,'partials',true,'leftSeed',sJo.Jo);

        
    elseif size(opt1.vRightSeed,1)~=0  %% forward with seed
        
        [o,oJv] = realizationOutput(x,u,v,sss,'partials',true,'vRightSeed',opt1.vRightSeed,'xRightSeed',opt1.xRightSeed,'uRightSeed',opt1.uRightSeed);
        
        J = oJv.J;
        J = bringVariablesMPI(J,sss.jobSchedule);
        
        if imMaster
        [s,Jac] = outputRisks(o,'eta',sss.eta,'partials',true,'oRightSeed',J);
        else
        Jac.J =nan;
        end
    else % full jacobian    % do it backards!  there are many more variables than outputs
        
        [o] = realizationOutput(x,u,v,sss,'partials',false);
        if imMaster       
        [s,sJo] = outputRisks(o,'eta',sss.eta,'partials',true);   
        else
        sJo.Jo =nan;    
        end
        [o,Jac] = realizationOutput(x,u,v,sss,'partials',true,'leftSeed',sJo.Jo);
        
    end
end

end


