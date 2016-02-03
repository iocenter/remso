function dh = calcDh(q,f,fref,nStages)
    dh  = pumpDhExplicit(q)*((f./fref).^2)*nStages;
end