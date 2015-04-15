function [o,oJac] = realizationOutput(x,u,v,outputFunction,varargin)
%  The full jacobian of this function is block diagonal w.r.t the
%  realizations.  --> Leave the jacobians in blocks


opt = struct('partials',false,'leftSeed',[],'vRightSeed',[],'xRightSeed',[],'uRightSeed',[]);
opt = merge_options(opt, varargin{:});

nR = numel(v);

u = repmat({u},nR,1);

rightSeedsGiven=true;
if ~(iscell(opt.vRightSeed) && size(opt.vRightSeed{1},1)~=0)
    rightSeedsGiven=false;
    opt.vRightSeed = cell(nR,1);
	opt.xRightSeed = cell(nR,1);
	opt.uRightSeed = cell(nR,1);
else
	opt.uRightSeed = repmat({opt.uRightSeed},nR,1);

end

if opt.partials && iscell(opt.leftSeed) && size(opt.leftSeed{1},2)~=0
    
    leftSeed = opt.leftSeed';
else
    leftSeed = cell(nR,1);    
end


if opt.partials
    [o,oJacD] = cellfun(@(of,xr,ur,vr,leftSeed,vRightSeed,xRightSeed,uRightSeed)of(xr,ur,vr,'gradients',true,'leftSeed',leftSeed,'vRightSeed',vRightSeed,'xRightSeed',xRightSeed,'uRightSeed',uRightSeed),...
                         outputFunction,x,u,v,leftSeed,opt.vRightSeed,opt.xRightSeed,opt.uRightSeed,'UniformOutput',false);
    if rightSeedsGiven
        oJac.J = cellfun(@(oJvr)oJvr.J,oJacD,'UniformOutput',false);
    else
        oJac.Jv = cellfun(@(oJvr)oJvr.Jv,oJacD,'UniformOutput',false)';
        oJac.Jx = cellfun(@(oJvr)oJvr.Jx,oJacD,'UniformOutput',false)';
        oJac.Ju = catAndSum(cellfun(@(oJvr)oJvr.Ju,oJacD,'UniformOutput',false));

    end
else
    [o] = cellfun(@(of,xr,ur,vr)of(xr,ur,vr),outputFunction,x,u,v,'UniformOutput',false);
    oJac = [];
end





end


function out = catAndSum(M)

dims = cellfun(@(x)size(x,2),M{1});
M = cellfun(@cell2mat,M,'UniformOutput',false);

if issparse(M{1})
    if isrow(M)
        M = M';
    end
    rows= size(M{1},1);
    blocks = numel(M);
    out = sparse( repmat(1:rows,1,blocks),1:rows*blocks,1)*cell2mat(M);
else
    out = sum(cat(3,M{:}),3);    
end

out = mat2cell(out,size(out,1),dims);

end