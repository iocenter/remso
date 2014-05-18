function varargout= simulateSystemResidual(x,u,ss,varargin)


opt = struct('gradients',false,'xLeftSeed',[],'vLeftSeed',[],'guessV',[],'xRightSeed',[],'uRightSeed',[]);
opt = merge_options(opt, varargin{:});

totalPredictionSteps = getTotalPredictionSteps(ss);
totalControlSteps = numel(u);

nx = numel(ss.state);
nv = ss.nv;
nu = numel(u{1});  %% assuming control dimension preserved

cs = cell(totalPredictionSteps,1);
v = cell(totalPredictionSteps,1);

if isempty(opt.guessV)
    opt.guessV = cell(totalPredictionSteps,1);
end

if opt.gradients
    JacStep = cell(1,totalPredictionSteps);
    if isempty(opt.xLeftSeed) && isempty(opt.xRightSeed)
        Jac.xJx = cell(totalPredictionSteps,totalPredictionSteps);
        Jac.xJu = cell(totalPredictionSteps,totalControlSteps);
        Jac.vJx = cell(totalPredictionSteps,totalPredictionSteps);
        Jac.vJu = cell(totalPredictionSteps,totalControlSteps);
    elseif ~isempty(opt.xRightSeed) && ~isempty(opt.xLeftSeed)
        error('not implemented')
    elseif ~isempty(opt.xRightSeed)
        Jac.xJ = cell(totalPredictionSteps,1);
        Jac.vJ = cell(totalPredictionSteps,1);
        
    elseif ~isempty(opt.xLeftSeed)
        Jac.Jx = cell(1,totalPredictionSteps,1);
        Jac.Ju = repmat({zeros(size(opt.xLeftSeed{1},1),nu)},1,totalControlSteps);
    end
else
    
end

%convergence = cell(1,totalPredictionSteps);   %%% can be output
converged = false(totalPredictionSteps,1);

%%% avoid "not sliced" variables, this additional variables do not represent copy of information!
ci = ss.ci;
if isempty(opt.xLeftSeed)
    xLeftSeed = cell(1,totalPredictionSteps);
else
    xLeftSeed = opt.xLeftSeed;
end
if isempty(opt.vLeftSeed)
    vLeftSeed = cell(1,totalPredictionSteps);
else
    vLeftSeed = opt.vLeftSeed;
end
if isempty(opt.xRightSeed)
    xRightSeed = cell(totalPredictionSteps,1);
else
    xRightSeed = opt.xRightSeed;
    %%% sort right hand sides!  %% Why .... because you have the off
    %%% diagonal term in the multipleshooting matrix!
    xRightSeed = [{zeros(nx,size(xRightSeed{1},2))};xRightSeed(1:end)];
end
uRightSeedSliced = cell(totalPredictionSteps,1);
if ~isempty(opt.uRightSeed)
    for k = 1:totalPredictionSteps
        uRightSeedSliced{k} = opt.uRightSeed{ss.ci(k)};    
    end	
end
usliced = cell(totalPredictionSteps,1);
for k = 1:totalPredictionSteps
        usliced{k} = u{ss.ci(k)};    
end	
step = ss.step;
gradientFlag = opt.gradients;
guessV = opt.guessV;
xStart = [ss.state;x(1:end-1)];

%diagonal terms
%parfor k = 1:totalPredictionSteps
for k = 1:totalPredictionSteps
    printCounter(1, totalPredictionSteps, k, 'ForwardSimMS')
    
    [cs{k},v{k},JacStep{k},convergence] = step{k}(xStart{k},usliced{k},...
        'gradients',gradientFlag,...
        'xLeftSeed',xLeftSeed{k},...
        'vLeftSeed',vLeftSeed{k},...
        'guessX',x{k},...
        'guessV',guessV{k},...
        'xRightSeeds',xRightSeed{k},...
        'uRightSeeds',uRightSeedSliced{k});
    
    
    converged(k) = convergence.converged;
    
    cs{k} = (cs{k}-x{k});
    
end
if opt.gradients
    if isempty(opt.xRightSeed) && isempty(opt.xLeftSeed)
        for k = 1:totalPredictionSteps
            [i] = feval(ci,k);
            Jac.xJu{k,i} = JacStep{k}.xJu;
            Jac.xJx{k,k} = -eye(nx);
            
            if nv >0
                Jac.vJu{k,i} = JacStep{k}.vJu;
            end
            %Jac.vJx{k,k} = -zeros(nv,nx);   it is assumed zero
            
            if k>1
                Jac.xJx{k,k-1} = JacStep{k}.xJx;
                if nv>0
                    Jac.vJx{k,k-1} = JacStep{k}.vJx;
                end
            elseif k == 1
            else
                error('what?')
            end
        end
        
    elseif ~isempty(opt.xRightSeed)
        for k = 1:totalPredictionSteps
            Jac.xJ{k} =  JacStep{k}.xJ - xRightSeed{k+1} ;
            if nv > 0
                Jac.vJ{k} = JacStep{k}.vJ;
            end
        end
        
    else %%% if ~isempty(opt.xLeftSeed)
        for k = 1:totalPredictionSteps
            [i] = ss.ci(k);
            Jac.Ju{1,i} = Jac.Ju{1,i}+JacStep{k}.Ju;
            Jac.Jx{1,k} = -xLeftSeed{k};
            if k>1
                Jac.Jx{1,k-1} = Jac.Jx{1,k-1}+JacStep{k}.Jx;
            elseif k == 1
            else
                error('what?')
            end
        end
    end
end

if ~all(converged)
    steps = 1:totalPredictionSteps;
    warning(strcat('System Residuals: Steps failed to converge:',num2str(steps(~converged))));
end

varargout{1} = cs;
varargout{2} = v;

if opt.gradients
    varargout{3} = Jac;
end
varargout{4} = converged;







end

