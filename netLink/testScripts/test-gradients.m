% close all
% p = rand(1000,1) * 100000;
%%%%%%%%%%%%%% 
%% Input %%%%%
%%%%%%%%%%%%%%
f = @superficialGasVelocity;
vars = {qgE, pV, zfactor, diameters, temperatures};
pert = [1e-6*meter^3/day; 1e-05*barsa; 1e-05; 1e-05; 1e-05];

%%%%%%%%%%%%%%%%%
%%% Algorithm %%%
%%%%%%%%%%%%%%%%%
for i=1:numel(vars)
    vars{i} = double(vars{i});   
end

fVal = f(vars{:});
gradAdi  = cell(1, numel(vars));  
gradPert = cell(1, numel(vars));  
for i=1:numel(vars)
    varsIt = vars;
    varsIt{i} = initVariablesADI(vars{i});
    tAdi = f(varsIt{:});
    
    
    grad = zeros(numel(fVal) , numel(double(vars{i})));
    
    for k = 1:numel(double(vars{i}))
        varsIt{i} = vars{i};        
        varsIt{i}(k) = vars{i}(k)+pert(i);
        
        grad(:,k) = (f(varsIt{:})-fVal)/pert(i);
    end
    gradAdi{i} = cell2mat(tAdi.jac);
    gradPert{i} = grad;    
end

a = full(cell2mat(gradAdi));
b = cell2mat(gradPert);
a-b


mrstVars = cell(1, numel(vars));
[mrstVars{:}] = initVariablesADI(vars{:});

fAdiFull = f(mrstVars{:});
c = full(cell2mat(fAdiFull.jac));

 
% plot(full((diag(tAdi.jac{:}))-grad)./grad)
