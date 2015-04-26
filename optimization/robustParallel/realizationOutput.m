function [o,oJac] = realizationOutput(x,u,v,sss,varargin)
%  The full jacobian of this function is block diagonal w.r.t the
%  realizations.  --> Leave the jacobians in blocks


opt = struct('partials',false,'leftSeed',[],'vRightSeed',[],'xRightSeed',[],'uRightSeed',[]);
opt = merge_options(opt, varargin{:});

jobSchedule = sss.jobSchedule;
ss = sss.ss;



gradients = opt.partials;
rightSeedsGiven=true;
vRightSeed = opt.vRightSeed;
xRightSeed = opt.xRightSeed;
uRightSeed = opt.uRightSeed;
if ~(iscell(uRightSeed) && size(uRightSeed{1},1)~=0)
    rightSeedsGiven=false;
    vRightSeed = createEmptyCompositeVar(jobSchedule);
	xRightSeed = createEmptyCompositeVar(jobSchedule);
end

leftSeed = opt.leftSeed;
leftSeedsGiven=true;
if isa(leftSeed,'Composite')
    % do no thing
elseif iscell(leftSeed)
    leftSeed = distributeVariables(leftSeed,jobSchedule);
else
    leftSeed = createEmptyCompositeVar(jobSchedule);
	leftSeedsGiven=false;
end


spmd
    [outputFunction] = getOutputF(ss);

    [o,J,Jx,Jv,Ju] = applyOutputFunction(outputFunction,x,u,v,leftSeed,vRightSeed,xRightSeed,uRightSeed,gradients,rightSeedsGiven);
    
    if gradients && leftSeedsGiven
        Ju = gop(@plus,Ju);
    end
end



o = bringVariables(o,jobSchedule);


oJac = [];
if gradients
    
    if rightSeedsGiven
        oJac.J = J;
    else
		Ju = Ju{1};

    	uDims = cellfun(@numel,u);
    	u = mat2cell(Ju,size(Ju,1),uDims);

        oJac.Jv = Jv;
        oJac.Jx = Jx;
        oJac.Ju = Ju;
        
    end

end


end


function out = catAndSum(M)
if isempty(M)
    out = 0;
else
    M = cellfun(@cell2mat,M,'UniformOutput',false);
    
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

end

function [o,J,Jx,Jv,Ju] = applyOutputFunction(outputFunction,x,u,v,leftSeed,vRightSeed,xRightSeed,uRightSeed,gradients,rightSeedsGiven)

J = [];
Jx = [];
Jv = [];
Ju = [];

[o,oJacD] = cellfun(@(...
    outputFunctionr,xr,vr,leftSeedr,vRightSeedr,xRightSeedr)...
    outputFunctionr(xr,u,vr,'gradients',gradients,'leftSeed',leftSeedr,'vRightSeeds',vRightSeedr,'xRightSeeds',xRightSeedr,'uRightSeeds',uRightSeed),...
    outputFunction, x, v ,leftSeed ,vRightSeed ,xRightSeed ,'UniformOutput',false);

if gradients
    if rightSeedsGiven
        J = cellfun(@(oJvr)oJvr.J,oJacD,'UniformOutput',false);
    else
        Jv = cellfun(@(oJvr)oJvr.Jv,oJacD,'UniformOutput',false);
        Jx = cellfun(@(oJvr)oJvr.Jx,oJacD,'UniformOutput',false);
        Ju = catAndSum(cellfun(@(oJvr)oJvr.Ju,oJacD,'UniformOutput',false));      
    end
end

end


function [outputF] = getOutputF(ss)
    outputF = cellfun(@(ssr)ssr.outputF,ss,'UniformOutput',false);
end


