function [o,oJv] = realizationOutput(v,outputFunction,varargin)
%  The full jacobian of this function is block diagonal w.r.t the
%  realizations.  --> Leave the jacobians in blocks


opt = struct('partials',false,'leftSeed',[],'vRightSeed',[]);
opt = merge_options(opt, varargin{:});

nR = numel(v);
rightSeedsGiven=true;
if ~(iscell(opt.vRightSeed) && size(opt.vRightSeed{1},1)~=0)
    rightSeedsGiven=false;
    opt.vRightSeed = cell(nR,1);
end

if opt.partials && iscell(opt.leftSeed) && size(opt.leftSeed{1},2)~=0
    
    leftSeed = opt.leftSeed';
else
    leftSeed = cell(nR,1);    
end


if opt.partials
    [o,oJvD] = cellfun(@(of,vr,leftSeed,vRightSeed)of(vr,'partials',true,'leftSeed',leftSeed,'vRightSeed',vRightSeed),outputFunction,v,leftSeed,opt.vRightSeed,'UniformOutput',false);
    if rightSeedsGiven
        oJv.J = cellfun(@(oJvr)oJvr.J,oJvD,'UniformOutput',false);
    else
        oJv.Jv = cellfun(@(oJvr)oJvr.Jv,oJvD,'UniformOutput',false)';
    end
else
    [o] = cellfun(@(of,vr)of(vr),outputFunction,v,'UniformOutput',false);
    oJv = [];
end





end

