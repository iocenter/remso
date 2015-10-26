

close all
p = rand(1000,1) * 100000;


pADI= initVariablesADI(p);

f = @(x)x.^(2);


pert = 1e-5;
grad = zeros(numel(double(p)),1);
for k = 1:numel(double(p))
    grad(k) = (f(p(k)+pert)-f(p(k)))/pert;
    
end


tAdi =   f(pADI);

plot(full((diag(tAdi.jac{1}))-grad)./grad)