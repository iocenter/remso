function [m,n] = linearCoefficients(x1,x2)
%% calculates linear coefficients passing through 2 points in R2
    m = (x2(:,2) - x1(:,2))./(x2(:,1) - x1(:,1));
    n = (x1(:,2)  -  m.*x1(:,1));    
end