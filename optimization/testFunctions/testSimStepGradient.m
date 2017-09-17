function [ maxError ] = testSimStepGradient(state,u,stepF,varargin)


opt = struct('pert',1e-6,'debug',false);
opt = merge_options(opt, varargin{:});



nx = numel(state);
nu = numel(u);


xRightSeeds = [eye(nx)     ,zeros(nx,nu)];
uRightSeeds = [zeros(nu,nx),eye(nu)];
[xp,v,Jac,convergence] = stepF(state,u,'gradients',true,'xRightSeeds',xRightSeeds,'uRightSeeds',uRightSeeds);

xJxFwd = Jac.xJ(:,1:nx); 
xJuFwd = Jac.xJ(:,nx+1:nx+nu);
nv = size(v,1);
if nv >0
    vJxFwd = Jac.vJ(:,1:nx); 
    vJuFwd = Jac.vJ(:,nx+1:nx+nu);
else
    vJxFwd = zeros(0,nx); 
    vJuFwd = zeros(0,nu);
    
end

xLeftSeed = [eye(nx);zeros(nv,nx)];
vLeftSeed = [zeros(nx,nv);eye(nv)];
[xp,v,Jac,convergence] = stepF(state,u,'gradients',true,'xLeftSeed',xLeftSeed,'vLeftSeed',vLeftSeed);
xJxAdj = Jac.Jx(1:nx,:); 
xJuAdj = Jac.Ju(1:nx,:);
nv = size(v,1);
if nv >0
    vJxAdj = Jac.Jx(nx+1:nx+nv,:); 
    vJuAdj = Jac.Ju(nx+1:nx+nv,:);
else
    vJxAdj = zeros(0,nx); 
    vJuAdj = zeros(0,nu);    
end

errorV = [norm(xJxFwd-xJxAdj),norm(xJuFwd-xJuAdj),norm(vJxFwd-vJxAdj),norm(vJuFwd-vJuAdj)];


JxP = zeros(size(xp,1),nx);
JuP = zeros(size(xp,1),nu);


vJxP = zeros(nv,nx);
vJuP = zeros(nv,nu);


for k = 1:nu
    
    uP = u;
    
    pertV = opt.pert;
    
    uP(k) = uP(k) + pertV;
    
    [xpP,vP] = stepF(state,uP);
    
    JuP(:,k) = (xpP-xp)/pertV;
    
    if nv >0
        vJuP(:,k) = (vP-v)/pertV;
    end
    if opt.debug
        if nv >0
            [k/nu,max(norm((xJuFwd(:,1:k)-JuP(:,1:k)),inf),norm((vJuFwd(:,1:k)-vJuP(:,1:k)),inf))]
        else
            [k/nu,max(norm((xJuFwd(:,1:k)-JuP(:,1:k)),inf))]
        end
    end
end


for k = 1:nx
    
    stateP = state;
    
    pertV = opt.pert;
    
    stateP(k) = stateP(k) + pertV;
    
    [xpP,vP] = stepF(stateP,u);
    
    
    JxP(:,k) = (xpP-xp)/pertV;
    
    if nv >0
        vJxP(:,k) = (vP-v)/pertV;
    end
    if opt.debug
        if nv >0
            [k/nx,max(norm((xJxFwd(:,1:k)-JxP(:,1:k)),inf),norm((vJxFwd(:,1:k)-vJxP(:,1:k)),inf))]
        else
            [k/nx,max(norm((xJxFwd(:,1:k)-JxP(:,1:k)),inf))]
        end
    end
end

if nv >0
    maxError = max([errorV,norm((xJxFwd-JxP),inf),norm((xJuFwd-JuP),inf),norm((vJxFwd-vJxP),inf),norm((vJuFwd-vJuP),inf)]);
else
    maxError = max([errorV,norm((xJxFwd-JxP),inf),norm((xJuFwd-JuP),inf)]);
    
end

end

