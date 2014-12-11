function  [ok,x]  = checkBounds( lbx,x,ubx,varargin)

opt = struct('chopp',false,'verbose',false);
opt = merge_options(opt, varargin{:});


ok = cellfun(@(l,v,u) all(l<=v) && all(v<=u),lbx,x,ubx);

if ~all(ok)
    if opt.chopp
        if opt.verbose
            warning('Variables out of bounds -> chopping varibles')
        end
        x = cellfun(@(l,v,u) min(max(l,v),u),lbx,x,ubx,'UniformOutput',false);
    else
        if opt.verbose
            warning('Variables out of bounds')
        end
    end
end


end

