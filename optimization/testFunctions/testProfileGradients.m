function [ ei,fi,vi,eFA ] = testProfileGradients(x,u,v,stepF,ci,state,varargin)


opt = struct('d',10,'pert',1e-5,'all',false);
opt = merge_options(opt, varargin{:});



nSteps = numel(stepF);


fi = zeros(nSteps,1);
ei = zeros(nSteps,1);
vi = zeros(nSteps,1);
eFA = zeros(nSteps,1);

for k = 1:nSteps
    
    cik = callArroba(ci,{k});
    
    if k ==1
        xStart = state;
    else
        xStart = x{k-1};
    end
    
    nx = numel(xStart);
    nu = numel(u{cik});
    n = nx+nu;
    
    if opt.all
        xRightSeed = sparse(nx,nu);
        uRightSeed = speye(nu);
        opt.d = nu;
    else
        rightSeed = [cell2mat(arrayfun(@(x)(orth(rand(n,opt.d))),1:floor(opt.d/n),'UniformOutput',false)),orth(rand(n,mod(opt.d,n)))] ;
        xRightSeed = rightSeed(1:nx,:);
        uRightSeed = rightSeed(nx+1:end,:);
    end
    
    
    [xs,vs,JacStep,~,simVars] = callArroba(stepF{k},{xStart,u{cik}},...
        'gradients',true,...
        'guessX',x{k},...
        'guessV',v{k},...
        'xRightSeeds',xRightSeed,...
        'uRightSeeds',uRightSeed);
    

    f = @(xx) callArroba(stepF{k},{xx(1:nx),xx(nx+1:nx+nu)},...
        'gradients',false,...
        'guessX',x{k},...
        'guessV',v{k});
    fm = @(xx) merge2Outs(f,xx);
    
	nx = numel(xs);
    nv = numel(vs);
    
    [g] = directionalGradFD(fm,[xStart;u{cik}],[xRightSeed;uRightSeed],opt.pert,nx+nv);
    
    g2 = [JacStep.xJ;JacStep.vJ];

    
    [maxV,i] = max(abs(g-g2));
    [maxV,j] = max(maxV);
    
    fi(k) = i(j);
    ei(k) = maxV;
    vi(k) = g2(i(j),j);
    
    

    if opt.all
        i = 1:(nx+nv);
        leftSeed = speye(nx+nv);
    else
        leftSeed = sparse(1:opt.d,i,1,opt.d,nx+nv);
    end
    xLeftSeed = leftSeed(:,1:nx);
    vLeftSeed = leftSeed(:,nx+1:nx+nv);

    [xs2,vs2,JacStepAdj] = callArroba(stepF{k},{xStart,u{cik}},...
        'gradients',true,...
        'guessX',x{k},...
        'guessV',v{k},...
        'xLeftSeed',xLeftSeed,...
        'vLeftSeed',vLeftSeed,...
        'simVars',simVars);
    
    assert(norm(xs-xs2,inf)<sqrt(eps))
    assert(norm(vs-vs2,inf)<sqrt(eps))
    
    F = g2(i,:);
    A = JacStepAdj.Jx*xRightSeed+JacStepAdj.Ju*uRightSeed;
    FA = F-A;
    eFA(k) = norm(FA,inf);

end


end


function [g] = directionalGradFD(f,x,dx,pert,ysize)

g = zeros(ysize,size(dx,2));

for k = 1:size(dx,2)
   g(:,k) = (f(x+dx(:,k)*pert)-f(x-dx(:,k)*pert))/(2*pert);
end


end


function [y] = merge2Outs(f,x)

[y1,y2] = f(x);

y = [y1;y2];

end
