function dp = calcDp(q, f, fref, nStages)        
    dh  = pumpDhExplicit(q)*((f./fref).^2)*nStages;
    mixDens = 900;
    dp  = (dh*norm(gravity)*mixDens)./barsa;

end
