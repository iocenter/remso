function [ varDistributed ] = distributeVariables( var,jobSchedule)
%
%  Given a cellarray variable defined for each time on the prediction
%  horizon, distributed the variables to workers according to jobSchedule


if isa(var,'Composite')
    varDistributed = var;
elseif isstruct(var)
    if isfield(var,'worker')
        varDistributed = var.worker;
    elseif isfield(var,'client')
        varDistributed = distributeToComposites(var.client,jobSchedule);
    else
        varDistributed = distributeToComposites(var,jobSchedule);        
    end
else  %% is a cell!
    varDistributed = distributeToComposites(var,jobSchedule);
end


end

function varDistributed = distributeToComposites(var,jobSchedule)
    
    varDistributed = Composite();
    for w = 1:length( varDistributed )
        varDistributed{w} = var(jobSchedule.work2Job{w});
    end
end

