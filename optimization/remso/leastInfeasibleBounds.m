function [lbs,ubs] = leastInfeasibleBounds(s,opt,withAlgs)
% modified bounds according to the least-infeasibilty bounds given by the QP
% algorithm
%
%    lbs = lb - s
%    ubs = ub + s
%

if isempty(s)
    ubs = {opt.ubx;
        opt.ubv};
    
    lbs = {opt.lbx;
        opt.lbv};
else
    if withAlgs
        ubs = {cellfun(@(x1,x2)(x1+x2),opt.ubx,s.x,'UniformOutput',false);
            cellfun(@(x1,x2)(x1+x2),opt.ubv,s.v,'UniformOutput',false)};
        
        lbs = {cellfun(@(x1,x2)(x1-x2),opt.lbx,s.x,'UniformOutput',false);
            cellfun(@(x1,x2)(x1-x2),opt.lbv,s.v,'UniformOutput',false)};
    else
        ubs = {cellfun(@(x1,x2)(x1+x2),opt.ubx,s.x,'UniformOutput',false)};
        
        lbs = {cellfun(@(x1,x2)(x1-x2),opt.lbx,s.x,'UniformOutput',false)};
        
    end
end
end