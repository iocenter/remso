
function plotDp(E, qo, qw, qg, p, incLiqFLow, incWcut, liqSteps, wcutSteps)  
    ql = qo+qw;

    qliqAxis = zeros(liqSteps, 1);
    wcutAxis = zeros(wcutSteps, 1);
    
    dp = zeros(liqSteps, wcutSteps, 1);
    for i=1:liqSteps
        qliqAxis(i) = ql +(i-1)*incLiqFLow;    
        for j=1:wcutSteps  
            if ql==0
                wcutAxis(j) = (j-1)*incWcut;
            else
                wcutAxis(j) = (qw/ql + (j-1)*incWcut);
            end
            
            dp(i,j) = dpBeggsBrill(E, qliqAxis(i)*(1-wcutAxis(j)), qliqAxis(i)*wcutAxis(j), qg, p);
        end
    end
    surf(qliqAxis/(meter^3/day), wcutAxis, dp/barsa);
    title('Beggs and Brill Correlation');
    xlabel('qliq (sm3/d)');
    ylabel('wcut (fraction)');
    zlabel('dp (bar)');

end