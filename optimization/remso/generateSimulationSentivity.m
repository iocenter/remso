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

if ~isempty(Jacs)
    m = arrayfun(@(JacI)size(JacI.Jx,1),Jacs);
    Jac = catJacs([Jacs;JacAct],xDims,vDims,uDims);
else
    m = 0;
    Jac=JacAct;
end

[~,Aact,~,~,~,~] = simulateSystemZ(u,x,v,ss,[],'simVars',simVars,'JacTar',Jac,'withAlgs',withAlgs);
Aact = cell2mat(Aact);
if m > 0
    mStart = cumsum([1;m(end-1)]);
    mEnd = cumsum(m);
    sensitivities = arrayfun(@(mS,mE)Aact(mS:mE,:),mStart,mEnd,'UniformOutput',false);
    sensitivities = [sensitivities;{Aact(mEnd(end)+1:end,:)}];
else
    sensitivities = {Aact(1:end,:)};
end
sensitivities = cellfun(@(s)mat2cell(s,size(s,1),uDims),sensitivities,'UniformOutput',false);

end

