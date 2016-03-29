function [ dpTotal ] = dpCVODES(Eout,  qo, qw, qg, p, varargin)
%dpCVODES calculates total pressure drop in a pipeline using CVODES

%% OBSERVATION:  There must NOT be ADI objects within Eout!

opt     = struct('dpFunction', @simpleDp, ...
    'backward', false,...
    'monitor',false,...
    'forwardGradient',false,...
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

% ---------------------
% CVODES initialization
% ---------------------
%% This is way too precise, consider relaxing in the actual application
InitialStep = 5;
options = CVodeSetOptions('UserData',data,...
    'RelTol',1.e-4,...
    'AbsTol',1.e-8,...
    'JacobianFn',@djacfn,...
    'MaxStep',integrationSteps,...
    'InitialStep',-InitialStep + 2*InitialStep*(opt.backward),...
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
if opt.backward
    xf =  pipeSizes;
else
    xf = -pipeSizes;
end


p0 = p;
CVodeInit(@rhsfn, 'BDF', 'Newton', x0, double(p0)/outputScaling, options);

if computeGradients && ~opt.forwardGradient
    CVodeAdjInit(150, 'Hermite');  %% TODO: study how to tune the magic number 150
end


% ------------------
% FSA initialization
% ------------------
if computeGradients && opt.forwardGradient
    
    Ns = 4;
    yS0 = [0,0,0,1];
    
    
    FSAoptions = CVodeSensSetOptions('method','Simultaneous','ErrControl', true); 
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
        'AbsTol',1.e-8,...
        'JacobianFn',@djacBfn);
    
    if opt.monitor
        mondataB = struct;
        mondataB.mode = 'both';
        optionsB = CVodeSetOptions(optionsB,...
            'MonitorFn','CVodeMonitorB',...
            'MonitorData', mondataB);
    end
    
    
    idxB = CVodeInitB(@rhsBfn, 'BDF', 'Newton', xf, 1, optionsB);
    
    optionsQB = CVodeQuadSetOptions('ErrControl',true,...
        'RelTol',1.e-4,'AbsTol',1.e-8);
    
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


dp = data.dpFunction(data.fixedParameters{:},data.variableParameters{:},p*data.outputScaling)/data.outputScaling;



flag = 0;
new_data = [];

end

% ===========================================================================

function [J, flag, new_data] = djacfn(x, p, fp, data)
% Dense Jacobian function


[p] = initVariablesADI(p);

dp = data.dpFunction(data.fixedParameters{:},data.variableParameters{:},p*data.outputScaling)/data.outputScaling;

%assert(fp == double(dp));

if isa(dp,'ADI')
    J = dp.jac{1};
else
    J = 0;  %% dp does not depend on p
end


flag = 0;
new_data = [];

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


[qo, qw, qg] = initVariablesADI(data.variableParameters{:});

dp = data.dpFunction(data.fixedParameters{:},qo, qw, qg ,p*data.outputScaling)/data.outputScaling;

qBd = -lambda'*[dp.jac{1}*data.paramScaling{1},...
                dp.jac{2}*data.paramScaling{2},...
                dp.jac{3}*data.paramScaling{3}];

flag = 0;
new_data = [];
end
% ===========================================================================

function [JB, flag, new_data] = djacBfn(t, y, yB, fyB, data)
% Backward problem Jacobian function

J = djacfn(t,y,fyB,data);
JB = -J';

flag = 0;
new_data = [];

end