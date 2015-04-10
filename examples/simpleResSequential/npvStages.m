function [o,oJv] = npvStages(v,varargin)
% we assume the (negative) of the cashflow (NPV) being accumulated in the last algebraic

opt = struct('partials',false,'leftSeed',[],'vRightSeed',[]);
opt = merge_options(opt, varargin{:});


vDims = cellfun(@numel,v);
o = cellfun(@(vi)vi(end),v);

oJv = [];
if opt.partials
    Jv = cellfun(@(vi,k)sparse(k,numel(vi),1,numel(v),numel(vi)),v',num2cell(1:numel(v)), 'UniformOutput',false);  %% this is the jacobian of the output function
    if size(opt.leftSeed,2)~=0
        
        oJv.Jv = mat2cell(opt.leftSeed*cell2mat(Jv),size(opt.leftSeed,1),vDims);
                
    elseif iscell(opt.vRightSeed) && size(opt.vRightSeed{1},1)~=0
             
        oJv.J = cell2mat(Jv)*cell2mat(opt.vRightSeed);
        
    else
        oJv.Jv = Jv;
    end
end


end

