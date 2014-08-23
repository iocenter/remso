function [ m,Jac,debugInfo ] = l1merit(f,dE,bE,ub,lb,rho,tau,varargin)
%  l1 merit function definition 
%
%  m = f + rho * norm(dE,1) + tau' * violation(lb,be,ub)
%
%  optional parameters:
%
%  gradients - true if the partial derivatives must be computed
%
%  fRightSeeds - For jacobain-vector product, f Right-hand-side
%  
%  dERightSeeds - For jacobain-vector product, dE Right-hand-side
%  
%  bERightSeeds - For jacobain-vector product, bE Right-hand-side
%
%  leftSeed = Left-hand-side for the vector-jacobian product.
%
%  Returns
%
%  m - merit value
%  
%  Jac - Jacobian
%
%  debugInfo - contain the value of the violations

opt = struct('gradients',false,'fRightSeeds',[],'dERightSeeds',[],'bERightSeeds',[],'leftSeed',1,'debug',false);
opt = merge_options(opt, varargin{:});


debugInfo = struct('f',0,'eq',0,'ineq',0);
debug = opt.debug;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


penalty = 0;
varN = numel(dE);

ubActC = cell(varN,1);
lbActiC = cell(varN,1);

for i = 1:varN
    
    dEi = dE{i};
    bEi = bE{i};
    ubi = ub{i};
    lbi = lb{i};
    taui = tau{i};
    spmd
        ubActi = cellfun(@gt,bEi,ubi,'UniformOutput',false);
        lbActi = cellfun(@lt,bEi,lbi,'UniformOutput',false);
        
        
        meqi  = rho * sum(cellfun(@sumAbs,dEi));
        mineqi = sum(cellfun(@boundPenaltyActiveUp,bEi,ubi,taui,ubActi))...
            +    sum(cellfun(@boundPenaltyActiveLow,bEi,lbi,taui,lbActi));
        
        p = gplus(meqi) + gplus(mineqi);
        
        
        if debug
            
            deqi =   max(cellfun(@maxAbs,dEi));
            dineqi = max( max(cellfun(@maxBoundUp,bEi,ubi,ubActi)),...
                max(cellfun(@maxBoundLow,bEi,lbi,lbActi)));
            
            deqi = gop(@max,deqi);
            dineqi = gop(@max,dineqi);
        end
        
        
    end
    ubActC{i} = ubActi;
    lbActiC{i} = lbActi;
    
    penalty = penalty + p{1};
    
    if debug
        debugInfo.eq = max(debugInfo.eq,deqi{1});
        debugInfo.ineq = max(debugInfo.ineq,dineqi{1});
    end
    
end

m = f + penalty;

if debug
    debugInfo.f = f;
end


if opt.gradients
    if ~(size(opt.fRightSeeds,1)==0)
        
        jp = 0;
        for i = 1:numel(dE)
            
            
            dEi = dE{i};
            dERightSeedsi = opt.dERightSeeds{i};
            taui = tau{i};
            ubActi = ubActC{i};
            lbActi = lbActiC{i};
            bERightSeedsi = opt.bERightSeeds{i};
            
            spmd                               
                jpC = rho * sum(cellfun(@eqLinePenaltyJac,dEi,dERightSeedsi))...
                    + sum(cellfun(@ineqLinePenaltyJac,taui,ubActi,lbActi,bERightSeedsi));
                
                jpC = gplus(jpC);
            end
            jp = jp + jpC{1};
        end
        Jac.J = opt.leftSeed*(opt.fRightSeeds + jp);       
    else
        error('not implemented')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end



function pg = eqLinePenaltyJac(di,rs)
pg = sign(di)'*rs;
end

function pg = ineqLinePenaltyJac(t,ua,la,rS)
pg = (t.*ua -t.*la)'*rS;
end

function ds = sumAbs(x)
ds = sum(abs(x));
end

function p = boundPenaltyActiveUp(v,b,t,a)
p = sum(t(a).*(v(a)-b(a)));
end

function p = boundPenaltyActiveLow(v,b,t,a)
p = sum(t(a).*(b(a)-v(a)));
end

function m = maxAbs(di)
m = max(abs(di));
end

function m = maxBoundUp(v,b,a)
m = max([0;(v(a)-b(a))]);
end

function m = maxBoundLow(v,b,a)
m = max([0;(b(a)-v(a))]);
end










