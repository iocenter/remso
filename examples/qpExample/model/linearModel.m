function [varargout] = linearModel(xStart,u,A,B,C,D,varargin)

opt = struct('gradients',false,'xLeftSeed',[],'vLeftSeed',[],'xRightSeeds',[],'guessV',[],'guessX',[],'uRightSeeds',[],'simVars',[],'W',0,'algFun',[]);
opt = merge_options(opt, varargin{:});

x = A*xStart + B * u + opt.W;
vi = C*xStart + D * u;


forwardStates =  x2states(x);
wellSols = v2wellSols(vi);

netSol = runNetwork(netSol, wellSol, forwardStates{step});



[no,JacAlg] = algFunc(forwardStates,wellSols,u);


v = [vi;no];






varargout{1} = x;
varargout{2} = v;
varargout{3} = [];


if opt.gradients

    xJx = A;
    xJu = B;
    vJx = C;
    vJu = D;
    
    if ~(size(opt.xLeftSeed,2)==0) && ~(size(opt.xRightSeeds,1)==0)
        error('not implemented')
    elseif (size(opt.xLeftSeed,2)==0) && (size(opt.xRightSeeds,1)==0)
        
        Jac.xJx = xJx;
        Jac.xJu = xJu;
        Jac.vJx = vJx;
        Jac.vJu = vJu;
        
    elseif ~(size(opt.xLeftSeed,2)==0)
        
        Jac.Jx = opt.xLeftSeed*xJx + opt.vLeftSeed * vJx;
        Jac.Ju = opt.xLeftSeed*xJu + opt.vLeftSeed * vJu;
        
    elseif ~(size(opt.xRightSeeds,1)==0)

        Jac.xJ = xJx*opt.xRightSeeds + xJu*opt.uRightSeeds;
        Jac.vJ = vJx*opt.xRightSeeds + vJu*opt.uRightSeeds;
        
    else
        error('what?')
    end
    varargout{3} = Jac;
end

varargout{4} = struct('converged',true); 
varargout{5} = [];    

