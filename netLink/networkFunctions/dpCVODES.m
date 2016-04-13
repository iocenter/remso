function [ dpTotal] = dpCVODES(Eout,  qo, qw, qg, p, varargin)
%dpCVODES calculates total pressure drop in a pipeline using CVODES

%% OBSERVATION:  There must NOT be ADI objects within Eout!

opt     = struct('dpFunction', @simpleDp, ...
    'monitor',false,...
    'forwardGradient',true,...
    'finiteDiff', false, ...
    'hasSurfaceGas', true, ...
    'pScale', 5*barsa,...
    'qlScale', 5*meter^3/day,...
    'qgScale', 100*(10*ft)^3/day);

opt     = merge_options(opt, varargin{:});

pipes = vertcat(Eout.pipeline);
pipeSizes = vertcat(pipes.len);
integrationSteps = vertcat(Eout.integrationStep);
numPipes = numel(Eout);

paramScaling = {opt.qlScale,opt.qlScale,opt.qgScale};
outputScaling = opt.pScale;
% and p is the initial condition and depends on these values


computeGradients = isa(qo,'ADI');  %% pass computePartials to this function to avoid this workaround.

data = struct();
data.dpFunction = opt.dpFunction;
data.hasSurfaceGas = opt.hasSurfaceGas;
data.qo = double(qo);
data.qg = double(qg);
data.qw = double(qw);
data.pipeSizes = pipeSizes;
data.Eout = Eout;
data.paramScaling = paramScaling;
data.outputScaling = outputScaling;
data.p = nan;
data.dp = nan;
data.numPipes = numPipes;


% ---------------------
% CVODES initialization
% ---------------------
%% This is way too precise, consider relaxing in the actual application

if opt.finiteDiff
    jacFun = @djacfd;    
else
    jacFun = @djacfn;    
end

InitialStep = min(integrationSteps./pipeSizes);
options = CVodeSetOptions('UserData',data,...
    'RelTol',1.e-4,...
    'AbsTol',1.e-2,...
    'JacobianFn',jacFun,...
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
xf = 1;
% xf = pipeSizes;


p0 = p;
CVodeInit(@rhsfn, 'Adams' , 'Functional', x0, double(p0)/outputScaling, options);

if computeGradients && ~opt.forwardGradient
    CVodeAdjInit(150, 'Hermite');  %% TODO: study how to tune the magic number 150
end


% ------------------
% FSA initialization
% ------------------
if computeGradients && opt.forwardGradient
    
    if opt.hasSurfaceGas
        yS0 = [zeros(numPipes,3*numPipes),eye(numPipes)];
    else
        yS0 = [zeros(numPipes,2*numPipes),eye(numPipes)];
    end
    
    
    
    if opt.finiteDiff
        sensFun = @rhsSfd;
    else
        sensFun = @rhsSfn;
    end
    
    FSAoptions = CVodeSensSetOptions('method','Staggered','ErrControl', false); 
    CVodeSensInit((3+opt.hasSurfaceGas)*numPipes, sensFun, yS0, FSAoptions);
    
    
    [status, x, pFinal, dpSens] = CVode(xf,'Normal');    
    
    if opt.hasSurfaceGas    
        jacPfinal = mat2cell((dpSens(:,           1:numPipes  )/paramScaling{1}*cell2mat(qo.jac) + ...
                              dpSens(:,  numPipes+1:numPipes*2)/paramScaling{2}*cell2mat(qw.jac) + ...
                              dpSens(:,2*numPipes+1:numPipes*3)/paramScaling{3}*cell2mat(qg.jac) + ...
                              dpSens(:,3*numPipes+1:numPipes*4)/outputScaling  *cell2mat(p.jac)...
                             )*outputScaling...
                             , numPipes, cellfun(@(x)size(x,2),qo.jac));
    else
        jacPfinal = mat2cell((dpSens(:,           1:numPipes  )/paramScaling{1}*cell2mat(qo.jac) + ...
                              dpSens(:,  numPipes+1:numPipes*2)/paramScaling{2}*cell2mat(qw.jac) + ...
                              dpSens(:,2*numPipes+1:numPipes*3)/outputScaling  *cell2mat(p.jac)...
                             )*outputScaling...
                             , numPipes, cellfun(@(x)size(x,2),qo.jac));
    end
    pFinal = ADI(pFinal*outputScaling,jacPfinal);
     
else
    [status, x, pFinal] = CVode(xf,'Normal');
    pFinal = pFinal*outputScaling;
end



if opt.monitor
    si = CVodeGetStats
end


if computeGradients && ~opt.forwardGradient
    
    if opt.finiteDiff
        jacFun = @djacBfd;
        quadFun = @quadBfd;
        rhsFun = @rhsBfd;
    else
        jacFun = @djacBfn;
        quadFun = @quadBfn;
        rhsFun = @rhsBfn;
    end
    
    %%TODO: check how to add 'SensDependent',true
    optionsB = CVodeSetOptions('UserData',data,...
            'RelTol',1.e-4,...
            'AbsTol',1.e-2,...
            'MaxStep',inf,...
            'InitialStep',-InitialStep,...%'SensDependent',true,...
        'JacobianFn',jacFun);
    
    if opt.monitor
        mondataB = struct;
        mondataB.mode = 'both';
        optionsB = CVodeSetOptions(optionsB,...
            'MonitorFn','CVodeMonitorB',...
            'MonitorData', mondataB);
    end
    
    
    idxB = CVodeInitB(rhsFun, 'Adams' , 'Functional', xf, ones(numPipes,1), optionsB);
    
    optionsQB = CVodeQuadSetOptions('ErrControl',false,...
        'RelTol',1.e-4,'AbsTol',1.e-2);
    
    if opt.hasSurfaceGas
        CVodeQuadInitB(idxB, quadFun, zeros(3*numPipes*numPipes,1), optionsQB);
    else
        CVodeQuadInitB(idxB, quadFun, zeros(2*numPipes*numPipes,1), optionsQB);        
    end
    % ----------------------------------------
    % Backward integration
    % ----------------------------------------
    
    
    [status,t,yB,dpSens] = CVodeB(0,'Normal');    

    dpSens = reshape(dpSens,numPipes,(2+opt.hasSurfaceGas)*numPipes);
    if opt.hasSurfaceGas
    
    jacPfinal = mat2cell((dpSens(:,1:numPipes)/paramScaling{1}*cell2mat(qo.jac) + ...
                          dpSens(:,numPipes+1:2*numPipes)/paramScaling{2}*cell2mat(qw.jac) + ...
                          dpSens(:,numPipes*2+1:3*numPipes)/paramScaling{3}*cell2mat(qg.jac) + ...
                         diag(yB)          /outputScaling  *cell2mat(p.jac)...
                         )*outputScaling...
                         ,numPipes, cellfun(@(x)size(x,2),qo.jac) );   
    else
    
    jacPfinal = mat2cell((dpSens(:,1)/paramScaling{1}*cell2mat(qo.jac) + ...
                          dpSens(:,2)/paramScaling{2}*cell2mat(qw.jac) + ...
                         yB          /outputScaling  *cell2mat(p.jac)...
                         )*outputScaling...
                         ,numPipes, cellfun(@(x)size(x,2),qo.jac) );           
    end
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
    dp = data.pipeSizes.*data.dpFunction(data.Eout,data.qo,data.qw,data.qg,p*data.outputScaling, data.hasSurfaceGas)/data.outputScaling;
    data.p=p;
    data.dp=dp;
    new_data = data;
end


flag = 0;

end

% ===========================================================================

function [J, flag, new_data] = djacfn(x, p, fp, data)
% Dense Jacobian function

if data.hasSurfaceGas
    [qo, qw, qg,p] = initVariablesADI(data.qo,data.qw,data.qg,p);
else
    [qo, qw,p] = initVariablesADI(data.qo,data.qw,p);
    qg = data.qg;
end

dp = data.pipeSizes.*data.dpFunction(data.Eout,qo, qw, qg ,p*data.outputScaling, data.hasSurfaceGas)/data.outputScaling;

%assert(fp == double(dp));

if isa(dp,'ADI')
    J = full(dp.jac{end});
else
    J = zeros(data.numPipes);  %% dp does not depend on p
end
data.p = double(p);
data.dp = dp;

flag = 0;
new_data = data;

end

% ===========================================================================
function [pSd, flag, new_data] = rhsSfn(t,p,fp,pS,data)
% Sensitivity right-hand side function


if data.hasSurfaceGas
    [qo, qw, qg,p] = initVariablesADI(data.qo,data.qw,data.qg,p);
else
    [qo, qw,p] = initVariablesADI(data.qo,data.qw,p);
    qg = data.qg;
end

dp = data.pipeSizes.*data.dpFunction(data.Eout,qo, qw, qg ,p*data.outputScaling, data.hasSurfaceGas)/data.outputScaling;                   

if data.hasSurfaceGas
pSd =   dp.jac{4}*pS+ [dp.jac{1}*data.paramScaling{1},...
                       dp.jac{2}*data.paramScaling{2},...
                       dp.jac{3}*data.paramScaling{3},...
                       zeros(data.numPipes)];
else
pSd =   dp.jac{4}*pS+ [dp.jac{1}*data.paramScaling{1},...
                       dp.jac{2}*data.paramScaling{2},...
                       zeros(data.numPipes)];
end
                   
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

if all(data.p == p) && isa(data.dp,'ADI')
    
    dp = data.dp;
    
else
    
    if data.hasSurfaceGas
        [qo, qw, qg] = initVariablesADI(data.qo,data.qw,data.qg);
    else
        [qo, qw] = initVariablesADI(data.qo,data.qw);
        qg = data.qg;
    end
    
    dp = data.pipeSizes.*data.dpFunction(data.Eout,qo, qw, qg ,p*data.outputScaling, data.hasSurfaceGas)/data.outputScaling;
end

if data.hasSurfaceGas
    
    qBd = -bsxfun(@times,lambda,[dp.jac{1}*data.paramScaling{1},...
                                 dp.jac{2}*data.paramScaling{2},...
                                 dp.jac{3}*data.paramScaling{3}]);
    
else
    
    qBd = -bsxfun(@times,lambda,[dp.jac{1}*data.paramScaling{1},...
                                 dp.jac{2}*data.paramScaling{2}]);
end
qBd = full(reshape(qBd,numel(qBd),1)); %%TODO: use a sparse matrix

            
flag = 0;
new_data = [];
end
% ===========================================================================

function [JB, flag, new_data] = djacBfn(t, y, yB, fyB, data)
% Backward problem Jacobian function

[J,flag,new_data] = djacfn(t,y,fyB,data);
JB = -J';

end


% ======================== FD-Forward ================================================

function [pSd, flag, new_data] = rhsSfd(t,p,fp,pS,data)
% Sensitivity right-hand side function with the finite differences method


[Jo, Jw, Jg, Jp] = dpGradFD(data.Eout, data.qo,data.qw,data.qg, p*data.outputScaling, data.hasSurfaceGas, [], [],  'dpFunction', data.dpFunction);  

Jp = bsxfun(@times,data.pipeSizes,Jp);
Jo = bsxfun(@times,data.pipeSizes,Jo)/data.outputScaling;
Jw = bsxfun(@times,data.pipeSizes,Jw)/data.outputScaling;
Jg = bsxfun(@times,data.pipeSizes,Jg)/data.outputScaling;

if data.hasSurfaceGas
pSd = Jp*pS + [Jo*data.paramScaling{1},...
               Jw*data.paramScaling{2},...
               Jg*data.paramScaling{3},...
               zeros(data.numPipes)];
else
pSd = Jp*pS + [Jo*data.paramScaling{1},...
               Jw*data.paramScaling{2},...
               zeros(data.numPipes)];    
end
pSd = full(pSd); %%TODO: use a sparse matrix

flag = 0;
new_data = [];

end

% ======================== FD-Backwards ================================================
function [yBd, flag, new_data] = rhsBfd(t, y, yB, data)
% Backward problem right-hand side function with finite diff
[JB, flag, new_data] = djacBfd(t, y, yB, [], data);

yBd = JB*yB;
end




function [qBd, flag, new_data] = quadBfd(x, p, lambda, data)
% Backward problem quadrature integrand function

if all(data.p == p) && isa(data.dp,'ADI')
    dp = data.dp;
    
    if data.hasSurfaceGas
        qBd = -lambda'*[dp.jac{1}*data.paramScaling{1},...
            dp.jac{2}*data.paramScaling{2},...
            dp.jac{3}*data.paramScaling{3}];
    else
        qBd = -lambda'*[dp.jac{1}*data.paramScaling{1},...
            dp.jac{2}*data.paramScaling{2},...
            dp.jac{3}*data.paramScaling{3}];
    end
else
    warning('Computing the Jp without need')

    [Jo, Jw, Jg, Jp] = dpGradFD(data.Eout, data.qo,data.qw,data.qg,p*data.outputScaling, data.hasSurfaceGas, [], [],  'dpFunction', data.dpFunction);   %%TODO: avoid extra computation of functions to obtain Jp

    Jo = bsxfun(@times,data.pipeSizes,Jo)/data.outputScaling;
    Jw = bsxfun(@times,data.pipeSizes,Jw)/data.outputScaling;
    
    if data.hasSurfaceGas
        Jg = bsxfun(@times,data.pipeSizes,Jg)/data.outputScaling;   
    
        qBd = -lambda'*[Jo*data.paramScaling{1},...
                        Jw*data.paramScaling{2},...
                        Jg*data.paramScaling{3}];
                    
    else    
        qBd = -lambda'*[Jo*data.paramScaling{1},...
                        Jw*data.paramScaling{2}];        
    end
    qBd = full(qBd);
end            
flag = 0;
new_data = [];
end
% ===========================================================================

function [JB, flag, new_data] = djacBfd(t, y, yB, fyB, data)
% Backward problem Jacobian function

[J,flag,new_data] = djacfd(t,y,fyB,data);
JB = -J';
end

function [J, flag, new_data] = djacfd(x, p, fp, data)
% Dense Jacobian function with the finite differences method

warning('Computing the Jo,Jw,Jg without need')


[Jo, Jw, Jg, Jp] = dpGradFD(data.Eout,data.qo,data.qw,data.qg, p*data.outputScaling, data.hasSurfaceGas, [], [],  'dpFunction', data.dpFunction);      %% TODO: avoid extra function evaluations to compute Jg, Jo, Jw  



J = bsxfun(@times,data.pipeSizes,Jp);
J = full(J);
flag = 0;
new_data = [];

end







