function [ sensitivities ] = generateSimulationSentivity(u,x,v,ss,simVars,Jacs,xDims,vDims,uDims,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

opt.eta = 0.9;

if numel(varargin) == 1
    activeSet = varargin{1};
else % second input option
	lowActive = varargin{1};
	upActive = varargin{2};
	activeSet.ub.x = upActive.x;
	activeSet.lb.x = lowActive.x;
	activeSet.ub.v = upActive.v;
	activeSet.lb.v = lowActive.v;
	activeSet.ub.s = upActive.s;
	activeSet.lb.s = lowActive.s;
end
withAlgs = true;


if ~isempty(Jacs)
    m = arrayfun(@(JacI)size(JacI.Js,1),Jacs);
else
    m = 0;
end

[lS,ns] = leftSeedSgen(activeSet,Jacs);
no = size(lS,1);
[~,Js] = realization2s(x,u,v,ss,'partials',true,'eta',opt.eta,'leftSeed',lS);


if sum(m) > 0
   Js = sumJacContribution(Js,Jacs,m,'Jx',no);
   Js = sumJacContribution(Js,Jacs,m,'Jv',no);
   Js = sumJacContributionS(Js,Jacs,m,'Ju',no);
    
end
Jsu = Js.Ju;

Js = rmfield(Js,'Ju');
%% Now Js contain the jacobians corresponding to the active set of s and Jacs


actCell = realizationActiveSetXV(activeSet);
[~,JacActXV,nlx,nux,nlv,nuv] = cellfun(@activeSet2TargetXV,actCell,'UniformOutput',false);


Js = realizationJacsXV(Js);

[JacActXVJs] = cellfun(@catJacsXV,JacActXV,Js',xDims,vDims,'UniformOutput',false);



[~,Aact,~,~,~,~] = cellfun(@(xr,vr,ssr,simVarsr,j)simulateSystemZ(u,xr,vr,ssr,[],'simVars',simVarsr,'JacTar',j,'withAlgs',withAlgs),x,v,ss,simVars,JacActXVJs,'UniformOutput',false);

%%% in Aact
% first  --> activeSet of X and V
% second --> activeSet of S
% third --> jacs

%how to return?
%lbx
%ubx
%lbv
%ubv
%lbs
%ubs
%
% Jacs

[alx,aux,alv,auv,sJ] = cellfun(@extractGradients,Aact,nlx,nux,nlv,nuv,'UniformOutput',false);



sJ = catAndSum(sJ)+cell2mat(Jsu);

Aact = [cell2mat([alx;aux;alv;auv]);sJ(1:ns,:)];


if m > 0
    mEnd = ns+cumsum(m);
    mStart = [ns+1;ns+1+mEnd(end-1)];
    sensitivities = arrayfun(@(mS,mE)sJ(mS:mE,:),mStart,mEnd,'UniformOutput',false);
    sensitivities = [sensitivities;{Aact}];
else
    sensitivities = {Aact};
end
sensitivities = cellfun(@(s)mat2cell(s,size(s,1),uDims),sensitivities,'UniformOutput',false);

end

function [lS,ns] = leftSeedSgen(act,Jacs)

actl = find(act.lb.s{1});
actu = find(act.ub.s{1});
nl = numel(actl);
nu = numel(actu);

ns = nl+nu;

lSact = sparse(1:ns,[actl;actu],[-ones(nl,1);ones(nu,1)],ns,numel(act.ub.s{1}));

if isempty(Jacs)
    lS = lSact;
else
    Js = vertcat(Jacs(:).Js);
    lS = [lSact;Js];
end

end

function Js = sumJacContribution(Js,Jacs,m,var,no)
mT = sum(m);
to = no-mT;
for k = 1:numel(m)
    from = to+1;
    to = to + m(k);
    if ~isempty(Jacs(k).(var))
        Js.(var) = sumJacContribS(Js.(var),Jacs(k).(var),from:to);
    end
end
assert(to==no)
end
function JM = sumJacContribS(JM,J,index)
f = @(JMr,Jr)sumJacContribSr(JMr,Jr,index);
JM = cellfun(f,JM,J,'UniformOutput',false);
end
function JM = sumJacContribSr(JMr,Jr,index)
JM = cellfun(@(JMrk,Jrk)subsAssSumM(JMrk,Jrk,index),JMr,Jr,'UniformOutput',false);
end
function JMrk = subsAssSumM(JMrk,Jrk,index)
JMrk(index,:) = JMrk(index,:)+Jrk;
end
function Js = sumJacContributionS(Js,Jacs,m,var,no)
mT = sum(m);
to = no-mT;
for k = 1:numel(m)
    from = to+1;
    to = to + m(k);
    if ~isempty(Jacs(k).(var))
        Js.(var) = sumJacContribSr(Js.(var),Jacs(k).(var),from:to);
    end
end
assert(to==no)
end

function actCell = realizationActiveSetXV(activeSet)
actCell = cellfun(@subsActiveSet,activeSet.lb.x,activeSet.ub.x,activeSet.lb.v,activeSet.ub.v,'UniformOutput',false);
end
function act = subsActiveSet(lbx,ubx,lbv,ubv)
act.lb.x = lbx;
act.ub.x = ubx;
act.lb.v = lbv;
act.ub.v = ubv;
end

function J = realizationJacsXV(Jac)
J = cellfun(@subsJacsXV,Jac.Jx,Jac.Jv,'UniformOutput',false);
end

function J = subsJacsXV(Jx,Jv)
J.Jv = Jv;
J.Jx = Jx;
end

function J = catJacsXV(JacActXV,Js,xDims,vDims)

J = catJacs([JacActXV,Js],xDims,vDims,[]);

end

function [alx,aux,alv,auv,sJ] = extractGradients(Aact,nlx,nux,nlv,nuv)
    Aact = cell2mat(Aact);
    
    index = [nlx;nux;nlv;nuv];
	final = cumsum(index);
    start = [1;1+final(1:end-1)];

	alx = Aact(start(1):final(1),:);
    aux = Aact(start(2):final(2),:);
 	alv = Aact(start(3):final(3),:);
    auv = Aact(start(4):final(4),:);   

    sJ = Aact(final(4)+1:end,:);

end


function out = catAndSum(M)

if any(cellfun(@issparse,M))
    if isrow(M)
        M = M';
    end
    rows= size(M{1},1);
    blocks = numel(M);
    out = sparse( repmat(1:rows,1,blocks),1:rows*blocks,1)*cell2mat(M);
else
    out = sum(cat(3,M{:}),3);    
end

end
