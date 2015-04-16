

mean = 0;
sigma = 1;
eta = 0.9;
tstar = norminv(eta,mean,sigma);




x = normrnd(mean,sigma,1000,1);


[cVar,cVarJAC,VaR] = cvar( x,eta );

VaR-tstar


%hist(x,numel(x)/10);

testJac = true;
if testJac
    
    minV = norminv(0.0001,mean,sigma);
    maxV = norminv(0.9999,mean,sigma);
    
    xt = minV:(maxV-minV)*0.001:maxV;
    cvarValue = zeros(1,numel(xt));
    cvarJac = zeros(1,numel(xt));
    
    for k=1:numel(xt)
        [cVari,cVarJACi,VaR] = cvar( [x;xt(k)],eta );
        cvarValue(k) = cVari;
        cvarJac(k) = cVarJACi(end);
    end
    
    figure
    p1 = subplot(3,1,1);
    plot(xt,cvarValue,'.-');
    p2 = subplot(3,1,2);
    plot(xt,cvarJac,'.-');
    p3 = subplot(3,1,3);
    hist(x,numel(x)/10);
    linkaxes([p1,p2,p3],'x')
    
    
    
    cvarJac2 = zeros(1,numel(x));
    
    pert = 0.0001;
    for k=1:numel(x)
        xP = x;
        xP(k) = xP(k) + pert;
        cvarJac2(k) = (cvar(xP,eta )- cVar)/pert ;
    end
    
    norm(cVarJAC-cvarJac2,inf)
end