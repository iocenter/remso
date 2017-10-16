function dp = calcDp(q, f, fref, nStages, varargin)        
    opt     = struct('mixDensity', 900);
    opt     = merge_options(opt, varargin{:});
    
    
    dh  = pumpDhExplicit(q).*((f./fref).^2).*nStages;    
    dp  = (dh.*norm(gravity).*opt.mixDensity)./barsa;
end
