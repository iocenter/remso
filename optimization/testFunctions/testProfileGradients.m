function [ ei,fi,vi ] = testProfileGradients(x,u,v,stepF,ci,state,varargin)


opt = struct('d',10,'pert',1e-5);
opt = merge_options(opt, varargin{:});


nSteps = numel(stepF);


fi = zeros(nSteps,1);
ei = zeros(nSteps,1);
vi = zeros(nSteps,1);

for k = 1:nSteps
    
    cik = callArroba(ci,{k});
    
    if k ==1
        xStart = state;
    else
        xStart = x{k-1};
    end
    
    nx = numel(xStart);
    nu = numel(u{cik});
    
    xRightSeed = rand(nx,opt.d)-0.5;
    uRightSeed = rand(nu,opt.d)-0.5;
    
    
    
    [xs,vs,JacStep] = callArroba(stepF{k},{xStart,u{cik}},...
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
    
    [g] = directionalGradFD(fm,[xStart;u{cik}],[xRightSeed;uRightSeed],opt.pert,numel(xs)+numel(vs));
    
    g2 = [JacStep.xJ;JacStep.vJ];
    
    [maxV,i] = max(abs(g-g2));
    [maxV,j] = max(maxV);
    
    fi(k) = i(j);
    ei(k) = maxV;
    vi(k) = g2(i(j),j);
    
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
