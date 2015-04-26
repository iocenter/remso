function  [ok,xN]  = checkBounds( lbx,x,ubx,varargin)

opt = struct('chopp',false,'verbose',false,'tol',1e-10);
opt = merge_options(opt, varargin{:});

cellBounds = iscell(lbx);

xN = x;
if cellBounds
    ok = cellfun(@(l,v,u) all(l<=v) && all(v<=u),lbx,x,ubx);
else
    if isnumeric(x) && isnumeric(ubx) 
        ok = all(lbx<=x) && all(x<=ubx);
    else
        ok = cellfun(@(v) all(lbx<=v) && all(v<=ubx),x);
    end
end


if ~all(ok)
    
    if cellBounds
        xR = cellfun(@(l,v,u) min(max(l,v),u),lbx,x,ubx,'UniformOutput',false);
    else
        if isnumeric(x) && isnumeric(ubx)
            xR = max(lbx,min(x,ubx));
        else
            xR = cellfun(@(v) min(max(lbx,v),ubx),x,'UniformOutput',false);
        end
    end
    if isnumeric(x)
        error = sum(abs(xR-x));
    else
        error = sum(cellfun(@(x1,x2)sum(abs(x1-x2)),xR,x));
    end
    
    if opt.chopp      
        xN = xR;     
        if opt.verbose
            if error > opt.tol
                warning(['Variables out of bounds -> chopping varibles, Violation norm1 = ',num2str(error)])
            end
        end       
    else
        if opt.verbose
            if error > opt.tol
                warning(['Variables out of bounds, Violation norm1 = ',num2str(error)])
            end
        end
    end
end


end

