function [ dpTotal ] = dpCVODES(Eout,  qo, qw, qg, p, varargin)
%dpCVODES calculates total pressure drop in a pipeline using CVODES

%% OBSERVATION:  There must NOT be ADI objects within Eout!

opt     = struct('dpFunction', @simpleDp, ...
    'monitor',false,...
    'forwardGradient',true,...
    'pScale', 5*barsa,...
    'qlScale', 5*meter^3/day,...
    'qgScale', 100*(10*ft)^3/day);

opt     = merge_options(opt, varargin{:});

pipes = vertcat(Eout.pipeline);
pipeSizes = vertcat(pipes.len);
integrationSteps = vertcat(Eout.integrationStep);


fixedParameters = {Eout};
variableParameters = {double(qo), double(qw), double(qg)};
paramScaling = {opt.qlScale,opt.qlScale,opt.qgScale};
outputScaling = opt.pScale;
% and p is the initial condition and depends on these values


computeGradients = isa(qo,'ADI');  %% pass computePartials to this function to avoid this workaround.

data = struct();
data.dpFunction = opt.dpFunction;
data.fixedParameters = fixedParameters;
data.variableParameters = variableParameters;
data.paramScaling = paramScaling;
data.outputScaling = outputScaling;
data.p = nan;
data.dp = nan;


% ---------------------
% CVODES initialization
% ---------------------
%% This is way too precise, consider relaxing in the actual application
InitialStep = integrationSteps;
options = CVodeSetOptions('UserData',data,...
    'RelTol',1.e-4,...
    'AbsTol',1.e-2,...
    'JacobianFn',@djacfn,...
    'MaxStep',inf,...
    'InitialStep',InitialStep,...
    'SensDependent',true);

if opt.monitor
    mondata = struct;
    mondata.mode = 'both';
    mondata.sol = true;
    mondata.sensi = true;
    options = CVodeSetOptions(options,'MonitorFn',@CVodeMonitor,'MonitorData',mondata);
end

%% OBS: the function is independent of the position in the pipeline
x0 = 0.0;
xf = pipeSizes;


p0 = p;
CVodeInit(@rhsfn, 'Adams' , 'Functional', x0, double(p0)/outputScaling, options);

if computeGradients && ~opt.forwardGradient
    CVodeAdjInit(150, 'Hermite');  %% TODO: study how to tune the magic number 150
end


% ------------------
% FSA initialization
% ------------------
if computeGradients && opt.forwardGradient
    
    Ns = 4;
    yS0 = [0,0,0,1];
    
    
    FSAoptions = CVodeSensSetOptions('method','Staggered','ErrControl', false); 
    CVodeSensInit(Ns, @rhsSfn, yS0, FSAoptions);
    
    
    [status, x, pFinal, dpSens] = CVode(xf,'Normal');
    
    jacPfinal = mat2cell((dpSens(1)/paramScaling{1}*cell2mat(qo.jac) + ...
                          dpSens(2)/paramScaling{2}*cell2mat(qw.jac) + ...
                          dpSens(3)/paramScaling{3}*cell2mat(qg.jac) + ...
                          dpSens(4)/outputScaling  *cell2mat(p.jac)...
                         )*outputScaling...
                         ,1, cellfun(@(x)size(x,2),qo.jac) );
    pFinal = ADI(pFinal*outputScaling,jacPfinal);
    
else
    [status, x, pFinal] = CVode(xf,'Normal');
    pFinal = pFinal*outputScaling;
end



if opt.monitor
    si = CVodeGetStats
end


if computeGradients && ~opt.forwardGradient
    
    %%TODO: check how to add 'SensDependent',true
    optionsB = CVodeSetOptions('UserData',data,...
            'RelTol',1.e-4,...
            'AbsTol',1.e-2,...
            'MaxStep',inf,...
            'InitialStep',-InitialStep,...%'SensDependent',true,...
        'JacobianFn',@djacBfn);
    
    if opt.monitor
        mondataB = struct;
        mondataB.mode = 'both';
        optionsB = CVodeSetOptions(optionsB,...
            'MonitorFn','CVodeMonitorB',...
            'MonitorData', mondataB);
    end
    
    
    idxB = CVodeInitB(@rhsBfn, 'Adams' , 'Functional', xf, 1, optionsB);
    
    optionsQB = CVodeQuadSetOptions('ErrControl',false,...
        'RelTol',1.e-4,'AbsTol',1.e-2);
    
    CVodeQuadInitB(idxB, @quadBfn, [0;0;0], optionsQB);
    
    % ----------------------------------------
    % Backward integration
    % ----------------------------------------
    
    
    [status,t,yB,dpSens] = CVodeB(0,'Normal');
    
    
    
    
    jacPfinal = mat2cell((dpSens(1)/paramScaling{1}*cell2mat(qo.jac) + ...
                          dpSens(2)/paramScaling{2}*cell2mat(qw.jac) + ...
                          dpSens(3)/paramScaling{3}*cell2mat(qg.jac) + ...
                         yB/outputScaling*cell2mat(p.jac)...
                         )*outputScaling...
                         ,1, cellfun(@(x)size(x,2),qo.jac) );
    pFinal = ADI(pFinal,jacPfinal);

end




% -----------
% Free memory
% -----------

CVodeFree;



dpTotal = pFinal-p0;




end




function [dp, flag, new_data] = rhsfn(x, p, data)

if p == data.p
    dp = double(data.dp);
    new_data = [];
else
    dp = data.dpFunction(data.fixedParameters{:},data.variableParameters{:},p*data.outputScaling)/data.outputScaling;
    data.p=p;
    data.dp=dp;
    new_data = data;
end


flag = 0;

end

% ===========================================================================

function [J, flag, new_data] = djacfn(x, p, fp, data)
% Dense Jacobian function


[qo, qw, qg,p] = initVariablesADI(data.variableParameters{:},p);

dp = data.dpFunction(data.fixedParameters{:},qo, qw, qg ,p*data.outputScaling)/data.outputScaling;

%assert(fp == double(dp));

if isa(dp,'ADI')
    J = dp.jac{4};
else
    J = 0;  %% dp does not depend on p
end
data.p = double(p);
data.dp = dp;

flag = 0;
new_data = data;

end



function [pSd, flag, new_data] = rhsSfn(t,p,fp,pS,data)
% Sensitivity right-hand side function


[qo, qw, qg,p] = initVariablesADI(data.variableParameters{:},p);

dp = data.dpFunction(data.fixedParameters{:},qo, qw, qg ,p*data.outputScaling)/data.outputScaling;


pSd = [dp.jac{4}]*pS+ [dp.jac{1}*data.paramScaling{1},...
                       dp.jac{2}*data.paramScaling{2},...
                       dp.jac{3}*data.paramScaling{3},0];



flag = 0;
new_data = [];

end

function [yBd, flag, new_data] = rhsBfn(t, y, yB, data)
% Backward problem right-hand side function


[JB, flag, new_data] = djacBfn(t, y, yB, [], data);

yBd = JB*yB;


end
% ===========================================================================

function [qBd, flag, new_data] = quadBfn(x, p, lambda, data)
% Backward problem quadrature integrand function

if data.p == p && isa(data.dp,'ADI')
    
    dp = data.dp;
    
else

    [qo, qw, qg] = initVariablesADI(data.variableParameters{:});

    dp = data.dpFunction(data.fixedParameters{:},qo, qw, qg ,p*data.outputScaling)/data.outputScaling;
end

qBd = -lambda'*[dp.jac{1}*data.paramScaling{1},...
                dp.jac{2}*data.paramScaling{2},...
                dp.jac{3}*data.paramScaling{3}];

            
flag = 0;
new_data = [];
end
% ===========================================================================

function [JB, flag, new_data] = djacBfn(t, y, yB, fyB, data)
% Backward problem Jacobian function

[J,flag,new_data] = djacfn(t,y,fyB,data);
JB = -J';

end