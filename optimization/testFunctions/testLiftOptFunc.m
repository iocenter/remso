function [ errorMax ] = testLiftOptFunc(x,u,ss,fT,withAlgs,varargin)

opt = struct('pert',1e-6,'testAdjoint',true,'testFwd',true);
opt = merge_options(opt, varargin{:});



[xp,v,Jft] = fT(x,u,ss,'gradients',true,'withAlgs',withAlgs);


xDims = cellfun(@(z)numel(z),x);
uDims = cellfun(@(z)numel(z),u);
vDims = cellfun(@(z)numel(z),v);

ADdxdx = cell2matFill(Jft.xJx,xDims,xDims');
ADdxdu = cell2matFill(Jft.xJu,xDims,uDims');
ADdvdx = cell2matFill(Jft.vJx,vDims,xDims');
ADdvdu = cell2matFill(Jft.vJu,vDims,uDims');



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
    xk = mat2cell(xk,xDims,1);
    [xpk,vk] = fT(xk,u,ss,'gradients',false);
    dxdx(:,k) = (cell2mat(xpk)-xp)/pertV;
    if withAlgs
        dvdx(:,k) = (cell2mat(vk)-v)/pertV;
    end
end
parfor k =1:numel(uV)
    %for k =1:numel(uV)
    uk = uV;
    
    uk(k) = uk(k) + pertV;
    uk = mat2cell(uk,uDims,1);
    [xpk,vk] = fT(x,uk,ss,'gradients',false);
    dxdu(:,k) = (cell2mat(xpk)-xp)/pertV;
    if withAlgs
        dvdu(:,k) = (cell2mat(vk)-v)/pertV;
    end
end


a = max(max(abs(dxdx-ADdxdx)));
b = max(max(abs(dxdu-ADdxdu)));
c = max(max(abs(dvdx-ADdvdx)));
d = max(max(abs(dvdu-ADdvdu)));

errorMax = max([a,b,c,d]);



if opt.testAdjoint

    
    adjVecX = rand(1,sum(xDims));
    adjVecV = rand(1,sum(vDims));
    
    adjVecX = mat2cell(adjVecX,1,xDims');
    if withAlgs
        adjVecV = mat2cell(adjVecV,1,vDims');
    else
        adjVecV = [];
    end
    [~,~,Jac] = fT(x,u,ss,'gradients',true,'xLeftSeed',adjVecX,'vLeftSeed',adjVecV);
    
    e = norm(cell2mat(Jac.Jx) - (cell2mat(adjVecX)*ADdxdx + cell2mat(adjVecV)*ADdvdx));
    f = norm(cell2mat(Jac.Ju) - (cell2mat(adjVecX)*ADdxdu + cell2mat(adjVecV)*ADdvdu));
    
    errorMax = max([errorMax,e,f]);
    
end

if opt.testFwd
    
    fwdVecX = rand(sum(xDims),1);
    fwdVecU = rand(sum(uDims),1);
    
    fwdVecX = mat2cell(fwdVecX,xDims,1);
    fwdVecU = mat2cell(fwdVecU,uDims,1);
    
    [~,~,Jac] = fT(x,u,ss,'gradients',true,'xRightSeed',fwdVecX,'uRightSeed',fwdVecU,'withAlgs',withAlgs);
    
    g = norm(cell2mat(Jac.xJ) - (ADdxdx*cell2mat(fwdVecX) + ADdxdu*cell2mat(fwdVecU)));
    if withAlgs
        h = norm(cell2mat(Jac.vJ) - (ADdvdx*cell2mat(fwdVecX) + ADdvdu*cell2mat(fwdVecU)));
    else
        h = 0;
    end
    
    errorMax = max([errorMax,g,h]);
end



end

