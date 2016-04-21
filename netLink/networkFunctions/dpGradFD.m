function [ Jo, Jw, Jg, Jp] = dpGradFD(E, qoE, qwE, qgE, pV,hasSurfaceGas,f_work,Z, varargin)
%dpGradFD calculates the pressure drop gradients using finite differences
%%TODO: be careful with flow rates close to zero.

opt     = struct('dpFunction', @simpleDp, ...
    'oilJac', true, ...
    'waterJac', true, ...
    'gasJac', true, ...
    'pressureJac', true);

opt     = merge_options(opt, varargin{:});

assert(opt.gasJac <= hasSurfaceGas); %% can only perturb gas when there is gas at the surface

if nargin < 6
    hasSurfaceGas = true;
end

if nargin < 7
    f_work = [];
end

if nargin < 8
    Z = [];
end
assert(numel(E) == numel(qoE)); %% number of edges should be equal to number of elements in oil flow vector

numPipes = numel(E);

pScale  = 5*barsa;
qlScale =  5*meter^3/day;
qgScale = 100*(10*ft)^3/day;

pert = 1e-4;
pertAxis = 2;

qlP = qlScale*pert;
qgP = qgScale*pert;
pP =  pScale*pert;

if opt.waterJac
    watPert = [qlP; -qlP];
else
    watPert = zeros(pertAxis,1);
end

if opt.oilJac
    oilPert = [qlP; -qlP];
else
    oilPert = zeros(pertAxis,1);
end
    
if opt.pressureJac    
    pPert   = [pP; -pP];
else
    pPert = zeros(pertAxis,1);
end

if opt.gasJac
    gasPert = [qgP; -qgP];
else
    gasPert = zeros(pertAxis,1);
end


if hasSurfaceGas      
    nPert = pertAxis*(opt.oilJac + opt.oilJac + opt.gasJac + opt.pressureJac);
    pertMatrix = repmat(blkdiag(gasPert, oilPert, watPert, pPert), numPipes, 1);
else
    nPert = pertAxis*(opt.oilJac + opt.oilJac + opt.pressureJac);
    pertMatrix = horzcat(zeros(numPipes*nPert,1), repmat(blkdiag(watPert, oilPert , pPert), numPipes, 1));    
end
pertMatrix( all(~pertMatrix,2), : ) = []; %% remove rows with only zeros

qFd =  reshape(repmat([qgE; qoE; qwE; pV], 1, nPert)', numPipes*nPert, 4) + pertMatrix;

if strcmp(func2str(opt.dpFunction),'dpBeggsBrillJDJ')
     [s,alpha,d,e,oil, rho_sc,s_in,s_out, T, T_in,T_out] = wrapperJDJ(E);
    
     TFd = reshape(repmat(T,1,nPert)',numPipes*nPert,1);
     alphaFd = reshape(repmat(alpha,1,nPert)', numPipes*nPert,1);
     dFd =  reshape(repmat(d,1,nPert)', numPipes*nPert,1);
     eFd = reshape(repmat(e,1,nPert)', numPipes*nPert,1);

     [dpVector, ~, ~] = Beggs_Brill_dpds(s,qFd(:,4),[],alphaFd,dFd,eFd,1,qFd(:,1),qFd(:,2),qFd(:,3),rho_sc,s-1,s+1,TFd,TFd, hasSurfaceGas, f_work, Z);
     
     dpPert = reshape(dpVector, nPert, numPipes)';
else   
    EFd = reshape(repmat(E,1,nPert)', numel(qoE)*nPert,1);
    dpPert = reshape(opt.dpFunction(EFd, qFd(:,2), qFd(:,3), qFd(:,1), qFd(:,4),'hasSurfaceGas', hasSurfaceGas ), nPert, numPipes)'; 
end

Jg = zeros(numPipes,1);
Jo = zeros(numPipes,1);
Jw = zeros(numPipes,1);
Jp = zeros(numPipes,1);

i = 1;

%% TODO: generalize this code for different number of pertubations        
if hasSurfaceGas    
    if opt.gasJac        
        Jg = (dpPert(:,i)-dpPert(:,i+1))/(2*qgP); 
        i = i+pertAxis;   
    end
    
    if opt.oilJac
        Jo = (dpPert(:,i)-dpPert(:,i+1))/(2*qlP);
        i = i + pertAxis;
    end
    
    
    if opt.waterJac        
        Jw = (dpPert(:,i)-dpPert(:,i+1))/(2*qlP);
        i = i + pertAxis;
    end
    
    if opt.pressureJac        
        Jp = (dpPert(:,i)-dpPert(:,i+1))/(2*pP);        
        i = i + pertAxis;
    end
else
    if opt.oilJac
        Jo = (dpPert(:,i)-dpPert(:,i+1))/(2*qlP);        
        i = i + pertAxis;
    end
    
    if opt.waterJac
        Jw = (dpPert(:,i)-dpPert(:,i+1))/(2*qlP);        
        i = i + pertAxis;
    end
    
    if opt.pressureJac
        Jp = (dpPert(:,i)-dpPert(:,i+1))/(2*pP);        
    end
end

index = 1:numPipes;
Jg = sparse(index,index,Jg);
Jo = sparse(index,index,Jo);
Jw = sparse(index,index,Jw);
Jp = sparse(index,index,Jp);


