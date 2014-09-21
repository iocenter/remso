function [ ] = writeSchedule(schedules,varargin)

opt = struct('fileName', 'schedule.inc');
opt = merge_options(opt, varargin{:});

fid = fopen(opt.fileName,'w');


for k = 1:numel(schedules)
    schedulek = schedules(k);
    ci = 1;
    currentControl = schedulek.step.control(ci);
    fprintControl(fid,schedulek.control(ci).W);
    
    fprintf(fid,'TSTEP\n%d ',schedulek.step.val(ci)/day);
    ci = ci+1;
    
    while ci <= numel(schedulek.step.control)
        
        if schedulek.step.control(ci) ~= currentControl
            fprintf(fid,'/\n\n');
            printControl(schedulek.control(ci));
            currentControl = schedulek.step.control(ci);
            fprintf(fid,'TSTEP\n');
        end
        fprintf(fid,'%d ',schedulek.step.val(ci)/day);
        
        ci = ci+1;
    end
    fprintf(fid,'/\n');
    
end

fclose(fid);

end

function [] = fprintControl(fid,W)

INJ = find(vertcat(W.sign) == 1);
PROD = find(vertcat(W.sign) == -1);


fprintf(fid,'WCONPROD\n');
for k = PROD'
    fprintf(fid,'%s %s %s %s 5* %d /\n',W(k).name,compi2Str(W(k).compi),status2Str(W(k).status),upper(W(k).type),scaledVal(W(k).val,W(k).type));
end
fprintf(fid,'/ \n');



fprintf(fid,'WCONINJE\n');
for k = INJ'
    fprintf(fid,'%s %s %s %s 2* %d /\n',W(k).name,compi2Str(W(k).compi),status2Str(W(k).status),upper(W(k).type),scaledVal(W(k).val,W(k).type));
end
fprintf(fid,'/ \n');



end

function [statusString] = status2Str(statusInt)
    if statusInt == 1
        statusString = 'OPEN';
    else
        statusString = 'SHUT';
    end

end


function [compiSTR] = compi2Str(compi)
    if     norm(compi - [1,0,0]) == 0
        compiSTR = 'WATER';
    elseif norm(compi - [0,1,0]) == 0
        compiSTR = 'OIL';
    elseif norm(compi - [0,0,1]) == 0
        compiSTR = 'GAS';
    elseif norm(compi - [1,1,0]) == 0
        compiSTR = 'LIQUID';
    else
        error('WHAT IS COMPI?')
    end

end

function sv = scaledVal(v,t)

switch upper(t)
    case 'BHP'
        sv = v/barsa;
    case 'RATE'
        sv = v*day;
    otherwise
        error('HOW SHOULD WE SCALE?')
end

end