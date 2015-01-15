function  [ok,xN]  = checkBounds( lbx,x,ubx,varargin)

opt = struct('chopp',false,'verbose',false,'tol',1e-10);
opt = merge_options(opt, varargin{:});

xN = x;
ok = cellfun(@(l,v,u) all(l<=v) && all(v<=u),lbx,x,ubx);

if ~all(ok)
    
    if opt.chopp

        xN = cellfun(@(l,v,u) min(max(l,v),u),lbx,x,ubx,'UniformOutput',false);
        

        if opt.verbose
            error = sum(cellfun(@(x1,x2)sum(abs(x1-x2)),xN,x));
            if error > opt.tol
                warning(['Variables out of bounds -> chopping varibles, Violation norm1 = ',num2str(error)])
            end
        end
        
        
    else
        if opt.verbose
            error = sum(cellfun(@(x1,x2)sum(abs(x1-x2)),xR,x));
            if error > opt.tol
                warning(['Variables out of bounds, Violation norm1 = ',num2str(error)])
            end
        end
    end
end


end

