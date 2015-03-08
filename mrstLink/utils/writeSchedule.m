function [ ] = writeSchedule(schedules,varargin)

opt = struct('fileName', 'schedule.inc','units','METRIC');
opt = merge_options(opt, varargin{:});

fid = fopen(opt.fileName,'w');


for k = 1:numel(schedules)
    schedulek = schedules(k);
    ci = 1;
    currentControl = schedulek.step.control(ci);
    fprintControl(fid,schedulek.control(ci).W,opt.units);
    
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

function [] = fprintControl(fid,W,units)

INJ = find(vertcat(W.sign) == 1);
PROD = find(vertcat(W.sign) == -1);


fprintf(fid,'WCONPROD\n');
for k = PROD'
    switch upper(W(k).type)
        case 'ORAT'
            fprintf(fid,'%s %s %s %d /\n',W(k).name,status2Str(W(k).status),upper(W(k).type),-scaledVal(W(k).val,W(k).type,units));
        case 'WRAT'
            fprintf(fid,'%s %s %s 1* %d /\n',W(k).name,status2Str(W(k).status),upper(W(k).type),-scaledVal(W(k).val,W(k).type,units));
        case 'GRAT'
            fprintf(fid,'%s %s %s 2* %d /\n',W(k).name,status2Str(W(k).status),upper(W(k).type),-scaledVal(W(k).val,W(k).type,units));
        case 'LRAT'
            fprintf(fid,'%s %s %s 3* %d /\n',W(k).name,status2Str(W(k).status),upper(W(k).type),-scaledVal(W(k).val,W(k).type,units));
        case 'CRAT'
            error('not implemented')
        case 'RESV'
            fprintf(fid,'%s %s %s 4* %d /\n',W(k).name,status2Str(W(k).status),upper(W(k).type),scaledVal(W(k).val,W(k).type,units));            
        case 'BHP'
            fprintf(fid,'%s %s %s 5* %d /\n',W(k).name,status2Str(W(k).status),upper(W(k).type),scaledVal(W(k).val,W(k).type,units));
        case 'THP'
            error('not implemented')
        otherwise
            error('not implemented')
    end
end
fprintf(fid,'/ \n');



fprintf(fid,'WCONINJE\n');
for k = INJ'
    compiString = compi2Str(W(k).compi);
    switch upper(W(k).type)
        case 'RATE'
            fprintf(fid,'%s %s %s %s %d /\n',W(k).name,compiString,status2Str(W(k).status),upper(W(k).type),scaledVal(W(k).val,W(k).type,units,compiString));
        case 'RESV'
            fprintf(fid,'%s %s %s %s 1* %d /\n',W(k).name,compiString,status2Str(W(k).status),upper(W(k).type),scaledVal(W(k).val,W(k).type,units,compiString));
        case 'BHP'
            fprintf(fid,'%s %s %s %s 2* %d /\n',W(k).name,compiString,status2Str(W(k).status),upper(W(k).type),scaledVal(W(k).val,W(k).type,units,compiString));
        case 'THP'
            error('not implemented')
        otherwise
            error('not implemented')
    end
    
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

function sv = scaledVal(v,t,units,compiString)

if nargin < 4
    compiString = 'LIQUID';
end

unitScale = unit_system(units);

switch upper(t)
    case 'BHP'
        sv = v/unitScale.press;
    case 'RATE'
        if strcmp(compiString,'GAS')
            sv = v/(unitScale.gasvol_s/unitScale.time);
        else  % must be any kind of liquid
            sv = v/(unitScale.liqvol_s/unitScale.time);            
        end
    case {'ORAT','WRAT','LRAT'}
        sv = v/(unitScale.liqvol_s/unitScale.time);
	case 'GRAT'
        sv = v/(unitScale.gasvol_s/unitScale.time);     
    otherwise
        error('HOW SHOULD WE SCALE?')
end

end

%% copied from mrst/modules/deckformat/deckinput/convertDeckUnits.m
function u = unit_system(rspec)
   metric = strcmp(rspec, 'METRIC');
   field  = strcmp(rspec, 'FIELD');
   lab    = strcmp(rspec, 'LAB');
   SI     = strcmp(rspec, 'SI');

   if sum([metric, field, lab, SI]) ~= 1,
      error(id('USys:Unknown'), ...
            'Input unit system must be either METRIC, FIELD, LAB, or SI.');
   end

   if metric,
      u = struct('length'   , meter             , ...
                 'time'     , day               , ...
                 'density'  , kilogram / meter^3, ...
                 'press'    , barsa             , ...
                 'concentr' , kilogram / meter^3, ... % Concentration
                 'compr'    , 1 / barsa         , ... % Compressibility
                 'viscosity', centi*poise       , ...
                 'perm'     , milli*darcy       , ...
                 'liqvol_s' , meter^3           , ... % Liquid vol, surf
                 'liqvol_r' , meter^3           , ... % Liquid vol, res
                 'gasvol_s' , meter^3           , ... % Gas vol, surf
                 'gasvol_r' , meter^3           , ... % Gas vol, res
                 'trans'    , centi*poise * meter^3 / (day * barsa));
   elseif field
      u = struct('length'   , ft                , ...
                 'time'     , day               , ...
                 'density'  , pound / ft^3      , ...
                 'press'    , psia              , ...
                 'concentr' , pound / stb       , ... % Concentration
                 'compr'    , 1 / psia          , ...
                 'viscosity', centi*poise       , ...
                 'perm'     , milli*darcy       , ...
                 'liqvol_s' , stb               , ...
                 'liqvol_r' , stb               , ...
                 'gasvol_s' , 1000 * ft^3       , ... % Mscf
                 'gasvol_r' , stb               , ...
                 'trans'    , centi*poise * stb / (day * psia));
   elseif lab,
      u = struct('length'   , centi*meter           , ...
                 'time'     , hour                  , ...
                 'density'  , gram / (centi*meter)^3, ...
                 'press'    , atm                   , ...
                 'concentr' , gram / (centi*meter)^3, ...
                 'compr'    , 1 / atm               , ...
                 'viscosity', centi*poise           , ...
                 'perm'     , milli*darcy           , ...
                 'liqvol_s' , (centi*meter)^3       , ...
                 'liqvol_r' , (centi*meter)^3       , ...
                 'gasvol_s' , (centi*meter)^3       , ...
                 'gasvol_r' , (centi*meter)^3       , ...
                 'trans'    , centi*poise * (centi*meter)^3 / (hour * atm));
   else
      % SI units.  MRST extension.  Idempotency.
      u = struct('length'   , 1, ...
                 'time'     , 1, ...
                 'density'  , 1, ...
                 'press'    , 1, ...
                 'concentr' , 1, ...
                 'compr'    , 1, ...
                 'viscosity', 1, ...
                 'perm'     , 1, ...
                 'liqvol_s' , 1, ...
                 'liqvol_r' , 1, ...
                 'gasvol_s' , 1, ...
                 'gasvol_r' , 1, ...
                 'trans'    , 1);
   end
end
