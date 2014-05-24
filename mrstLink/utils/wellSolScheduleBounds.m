function [wellSolMax,wellSolMin] = wellSolScheduleBounds(wellSol,maxProd,maxInj,minProd,minInj)
%
% fill a wellSol mock object with maximums and minimums according to the
% input parameters
%

wellSolMax = wellSol;
wellSolMin = wellSol;

for w = 1:numel(wellSol)
    switch wellSol(w).sign
        case -1
            wellSolMax(w).bhp = maxProd.BHP;
            wellSolMax(w).qOs = -minProd.ORAT;
            wellSolMax(w).qWs = -minProd.WRAT;
            wellSolMin(w).bhp = minProd.BHP;
            wellSolMin(w).qOs = -maxProd.ORAT;
            wellSolMin(w).qWs = -maxProd.WRAT;
        case 1
            wellSolMax(w).bhp = maxInj.BHP;
            wellSolMax(w).qOs = maxInj.RATE;
            wellSolMax(w).qWs = maxInj.RATE;
            wellSolMin(w).bhp = minInj.BHP;
            wellSolMin(w).qOs = minInj.RATE;
            wellSolMin(w).qWs = minInj.RATE;
        otherwise
            error('what');
    end
    
end




end