function [ vals, type, control ] = getValuesFromSchedule(schedule)   %%% or schedule2CellControls
%
% extract schedule control information 
%

vals = zeroVals(schedule);
control = zeroStrings(schedule);
type = zeroStrings(schedule);

for k = 1:numel(schedule.control)
    
    inj    = schedule.control(k).WCONINJE;
    numInj = size(inj, 1);
    
    for ki = 1:numInj
        [vals{k}(ki),control{k}{ki}] = getInjectorValue(inj(ki,:));
         type{k}{ki} = 'inj';
    end
    
    prod    = schedule.control(k).WCONPROD;
    numProd = size(prod,1);
    
    for kp = 1:numProd
        [vals{k}(numInj+kp),control{k}{numInj+kp}] = getProducerValue(prod(kp,:));
         type{k}{numInj+kp} = 'prod';
    end
end

end

%--------------------------------------------------------------------------

function [val,cntrMode] = getInjectorValue(inj)
cntrMode = inj{4};
switch cntrMode
    case 'RATE'
        inx = 5;
    case 'RESV'
        inx = 6;
    case 'BHP'
        inx = 7;
    otherwise
        error(['Cannot handle control mode: ', cntrMode]);
end
val = inj{inx};
end

function [val,cntrMode] = getProducerValue(prod)
cntrMode = prod{3};
switch cntrMode
    case 'ORAT'
        inx = 4;
    case 'WRAT'
        inx = 5;
    case 'LRAT'
        inx = 6;
    case 'RESV'
        inx = 8;
    case 'BHP'
        inx = 9;
    otherwise
        error(['Cannot handle control mode: ', cntrMode]);
end

val = prod{inx};

end

function zz = zeroVals(schedule)
zz = cell(1, numel(schedule.control));
for k = 1:numel(zz)
    zz{k} = zeros( size(schedule.control(k).WCONINJE, 1) + ...
        size(schedule.control(k).WCONPROD, 1) , 1);
end

end

function zz = zeroStrings(schedule)
zz = cell(1, numel(schedule.control));
for k = 1:numel(zz)
    zz{k} = cell( size(schedule.control(k).WCONINJE, 1) + ...
        size(schedule.control(k).WCONPROD, 1) , 1);
end

end

