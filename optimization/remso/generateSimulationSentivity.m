function [ sensitivities ] = generateSimulationSentivity(u,x,v,ss,simVars,Jacs,xDims,vDims,uDims,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if numel(varargin) == 1
    activeSet = varargin{1};
else % second input option
	lowActive = varargin{1};
	upActive = varargin{2};
	activeSet.ub.x = upActive.x;
	activeSet.lb.x = lowActive.x;
	activeSet.ub.v = upActive.v;
	activeSet.lb.v = lowActive.v;
end
withAlgs = true;


[~,JacAct ] = activeSet2TargetXV(activeSet);

Jacsxv = rmfield(Jacs,'Ju');

if ~isempty(Jacsxv)
    m = arrayfun(@(JacI)size(JacI.Jx,1),Jacsxv);
    Jac = catJacs([Jacsxv;JacAct],xDims,vDims,uDims);
else
    m = 0;
    Jac=JacAct;
end

[~,Aact] = simulateSystemZ(u,x,v,ss,[],'simVars',simVars,'JacTar',Jac,'withAlgs',withAlgs);
Aact = cell2mat(Aact);
if sum(m) > 0
    mStart = cumsum([1;m(1:end-1)]);
    mEnd = cumsum(m);
    mi = (1:numel(m))';
    sensitivities = arrayfun(@(mS,mE,mii)Aact(mS:mE,:)+cell2mat(Jacs(mii).Ju),mStart,mEnd,mi,'UniformOutput',false);
    sensitivities = [sensitivities;{Aact(mEnd(end)+1:end,:)}];
else
    sensitivities = {Aact};
end
sensitivities = cellfun(@(s)mat2cell(s,size(s,1),uDims),sensitivities,'UniformOutput',false);

end

