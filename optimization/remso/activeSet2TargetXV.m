function [ t,Jac ] = activeSet2TargetXV(uDims,withAlgs,varargin )

if nargin == 3
    activeSet = varargin{1};
else % second input option
	lowActive = varargin{1};
	upActive = varargin{2};
	activeSet.ub.x = upActive.x;
	activeSet.lb.x = lowActive.x;
	if withAlgs
        activeSet.ub.v = upActive.v;
        activeSet.lb.v = lowActive.v;
	end 
end
    
    
xDims = cellfun(@numel,activeSet.lb.x);
vDims = cellfun(@numel,activeSet.lb.v);


Jxjm = find(cell2mat(activeSet.lb.x));
Jxjp = find(cell2mat(activeSet.ub.x));

Jvjm = find(cell2mat(activeSet.lb.v));
Jvjp = find(cell2mat(activeSet.ub.v));

outsxm = numel(Jxjm);
outsxp = numel(Jxjp);

outsvm = numel(Jvjm);
outsvp = numel(Jvjp);


outsX = outsxm+outsxp;
outsV = outsvm+outsvp;

outs = (outsX+outsV);


Jx = mat2cell(sparse(1:outsX,[Jxjm;Jxjp],[-ones(outsxm,1);ones(outsxp,1)],outs,sum(xDims)),outs,xDims);
Jv = mat2cell(sparse(outsX+(1:outsV),[Jvjm;Jvjp],[-ones(outsvm,1);ones(outsvp,1)],outs,sum(vDims)),outs,vDims);
Ju = mat2cell(sparse(outs,sum(uDims)),outs,uDims);


t = @(x,u,v) target(x,u,v,activeSet,Jx,Jv,Ju);

Jac.Jx = Jx;
Jac.Jv = Jv;
Jac.Ju = Ju;

end



function [obj,Jac] = target(x,v,activeSet,Jx,Jv,Ju)
% TODO: for generality, implement the treatment of the options

obj = cell2mat([cellfun(@(xi,acti)-x(act),x,activeSet.lb.x,'UniformOutput',false);
                cellfun(@(xi,acti) x(act),x,activeSet.ub.x,'UniformOutput',false); 
                cellfun(@(xi,acti)-x(act),v,activeSet.lb.v,'UniformOutput',false);
                cellfun(@(xi,acti) x(act),v,activeSet.ub.v,'UniformOutput',false)]);   

                       
Jac.Jx = Jx;
Jac.Jv = Jv;
Jac.Ju = Ju;
        
            
end
