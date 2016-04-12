function [ Jo, Jw, Jg, Jp] = dpGradFD(E, qoE, qwE, qgE, pV,hasSurfaceGas,f_work,Z, varargin)
%dpGradFD calculates the pressure drop gradients using finite differences
%%TODO: be careful with flow rates close to zero.

 opt     = struct('dpFunction', @simpleDp);
 opt     = merge_options(opt, varargin{:});
 
if nargin < 6
    hasSurfaceGas = true;
end

if nargin < 7
    f_work = [];
end

if nargin < 8
   Z = [];
end
 
pScale  = 5*barsa;
qlScale =  5*meter^3/day;
qgScale = 100*(10*ft)^3/day;

pert = 1e-4;

qlP = qlScale*pert;
qgP = qgScale*pert;
pP =  pScale*pert;


if hasSurfaceGas     
    nPert = 8;
    pertMatrix = repmat(blkdiag([qgP; -qgP], [qlP; -qlP],[qlP; -qlP], [pP; -pP]), numel(qoE), 1);
else
    nPert = 6;
    pertMatrix = horzcat(zeros(numel(qoE)*nPert,1), repmat(blkdiag([qlP; -qlP],[qlP; -qlP], [pP; -pP]), numel(qoE), 1));    
end

qFd =  reshape(repmat([qgE; qoE; qwE; pV], 1, nPert)', numel(qoE)*nPert, 4) + pertMatrix;

if strcmp(func2str(opt.dpFunction),'dpBeggsBrillJDJ')
     [s,alpha,d,e,oil, rho_sc,s_in,s_out, T, T_in,T_out] = wrapperJDJ(E);
    
     TFd = reshape(repmat(T,1,nPert)', numel(qoE)*nPert,1);
     alphaFd = reshape(repmat(alpha,1,nPert)', numel(qoE)*nPert,1);
     dFd =  reshape(repmat(d,1,nPert)', numel(qoE)*nPert,1);
     eFd = reshape(repmat(e,1,nPert)', numel(qoE)*nPert,1);

     [dpVector, ~, ~] = Beggs_Brill_dpds(s,qFd(:,4),[],alphaFd,dFd,eFd,1,qFd(:,1),qFd(:,2),qFd(:,3),rho_sc,s-1,s+1,TFd,TFd, hasSurfaceGas, f_work, Z);
     
     dpPert = reshape(dpVector, nPert, numel(qoE))';
else   
    %% TODO: include optional parameters such as hasSurfaceGas
    %% TODO: Check if passing only E without extra copies shows an desired behavior.   
    EFd = reshape(repmat(E,1,nPert)', numel(qoE)*nPert,1);
    
    dpPert = reshape(opt.dpFunction(EFd, qFd(:,2), qFd(:,3), qFd(:,1), qFd(:,4),'hasSurfaceGas', hasSurfaceGas ), nPert, numel(qoE))'; 
end

if hasSurfaceGas
    Jg = (dpPert(:,1)-dpPert(:,2))/(2*qgP);
    Jo = (dpPert(:,3)-dpPert(:,4))/(2*qlP);
    Jw = (dpPert(:,5)-dpPert(:,6))/(2*qlP);
    Jp = (dpPert(:,7)-dpPert(:,8))/(2*pP);
else
    Jg = zeros(numel(qoE),1);
    Jo = (dpPert(:,1)-dpPert(:,2))/(2*qlP);
    Jw = (dpPert(:,3)-dpPert(:,4))/(2*qlP);
    Jp = (dpPert(:,5)-dpPert(:,6))/(2*pP);
end
