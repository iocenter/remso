function [ errorMax ] = testLiftOptFunc(x,u,ss,fT,varargin)

opt = struct('pert',1e-6,'testAdjoint',true,'testFwd',true);
opt = merge_options(opt, varargin{:});



[xp,v,Jft] = fT(x,u,ss,'gradients',true);


nx = numel(x{1});
nu = numel(u{1});
nv = numel(v{1});

ADdxdx = cell2matFill(Jft.xJx,[nx,nx]);
ADdxdu = cell2matFill(Jft.xJu,[nx,nu]);
ADdvdx = cell2matFill(Jft.vJx,[nv,nx]);
ADdvdu = cell2matFill(Jft.vJu,[nv,nu]);



uV = cell2mat(u);
xV = cell2mat(x);
xp = cell2mat(xp);
v = cell2mat(v);

dxdx = zeros(numel(xp),numel(xV));
dxdu = zeros(numel(xp),numel(uV));

dvdx = zeros(numel(v),numel(xV));
dvdu = zeros(numel(v),numel(uV));

pertV = opt.pert;

parfor k =1:numel(xV)
    %for k =1:numel(xV)
    
    xk = xV;
    
    
    xk(k) = xk(k) + pertV;
    xk = toStructuredCells( xk,nx);
    [xpk,vk] = fT(xk,u,ss,'gradients',false);
    dxdx(:,k) = (cell2mat(xpk)-xp)/pertV;
    if nv >0
        dvdx(:,k) = (cell2mat(vk)-v)/pertV;
    end
end
parfor k =1:numel(uV)
    %for k =1:numel(uV)
    uk = uV;
    
    uk(k) = uk(k) + pertV;
    uk = toStructuredCells( uk,nu);
    [xpk,vk] = fT(x,uk,ss,'gradients',false);
    dxdu(:,k) = (cell2mat(xpk)-xp)/pertV;
    if nv >0
        dvdu(:,k) = (cell2mat(vk)-v)/pertV;
    end
end


a = norm( dxdx-ADdxdx);
b = norm( dxdu-ADdxdu);
c = norm( dvdx-ADdvdx);
d = norm( dvdu-ADdvdu);

errorMax = max([a,b,c,d]);



if opt.testAdjoint
    totalPredictedSteps = numel(x);
    
    adjVecX = rand(1,nx*totalPredictedSteps);
    adjVecV = rand(1,nv*totalPredictedSteps);
    
    adjVecX = toStructuredCells(adjVecX,nx,'T',true);
    if nv>0
        adjVecV = toStructuredCells(adjVecV,nv,'T',true);
    else
        adjVecV = repmat({zeros(1,0)},1,totalPredictedSteps);
    end
    [~,~,Jac] = fT(x,u,ss,'gradients',true,'xLeftSeed',adjVecX,'vLeftSeed',adjVecV);
    
    e = norm(cell2mat(Jac.Jx) - (cell2mat(adjVecX)*ADdxdx + cell2mat(adjVecV)*ADdvdx));
    f = norm(cell2mat(Jac.Ju) - (cell2mat(adjVecX)*ADdxdu + cell2mat(adjVecV)*ADdvdu));
    
    errorMax = max([errorMax,e,f]);
    
end

if opt.testFwd
    
    totalPredictedSteps = numel(x);
    totalControlSteps = numel(x);
    
    fwdVecX = rand(nx*totalPredictedSteps,1);
    fwdVecU = rand(nu*totalControlSteps,1);
    
    fwdVecX = toStructuredCells(fwdVecX,nx);
    fwdVecU = toStructuredCells(fwdVecU,nu);
    
    [~,~,Jac] = fT(x,u,ss,'gradients',true,'xRightSeed',fwdVecX,'uRightSeed',fwdVecU);
    
    g = norm(cell2mat(Jac.xJ) - (ADdxdx*cell2mat(fwdVecX) + ADdxdu*cell2mat(fwdVecU)));
    if nv>0
        h = norm(cell2mat(Jac.vJ) - (ADdvdx*cell2mat(fwdVecX) + ADdvdu*cell2mat(fwdVecU)));
    else
        h = 0;
    end
    
    errorMax = max([errorMax,g,h]);
end



end

