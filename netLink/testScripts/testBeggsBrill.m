%clear all
%close all

mrstModule add deckformat
mrstModule add ad-fi

format long g

n= 100;

s= rand(1,1);
alpha = -pi/2 + 100*pi*rand(n,1);
d = 1*inch + 20*inch*rand(n,1);
e = 1e-5 + 3e-5*rand(n,1);
T = 0 + 125*rand(n,1);


R_sbBounds = [0,254];
qoBounds = [0;1000*(meter^3/day)];
qwBounds = [0;1000*(meter^3/day)];
qgBounds = [0;10000*(meter^3/day)];  
pBounds = [2*barsa;500*barsa];  

q_o_sc = qoBounds(1) + (qoBounds(2)-qoBounds(1))*rand(n,1);
q_w_sc = qwBounds(1) + (qwBounds(2)-qwBounds(1))*rand(n,1);
q_g_sc = (q_o_sc*R_sbBounds(1) + q_o_sc.*(R_sbBounds(2)-R_sbBounds(1)).*rand(n,1));
pV =  pBounds(1) + ( pBounds(2)- pBounds(1))*rand(n,1);

rho_sc = [0.8000;906.6450;1000];

hasSurfaceGas = true;


addpath(genpath('/home/codas/Desktop/Transfer/Fluid properties'));
addpath(genpath('/home/codas/Desktop/Transfer/Pipe flow'));
addpath(genpath('/home/codas/Desktop/Transfer/Conversion factors'));


tic
dpVector = zeros(n,1);
for k =1:n
    [dp] = Beggs_Brill_dpds(s,pV(k),[],alpha(k),d(k),e(k),1,[q_g_sc(k),q_o_sc(k),q_w_sc(k)],rho_sc',s-1,s+1,T(k),T(k));
    dpVector(k) = dp(1);
end
forTime = toc



addpath(genpath('../../netLink/fluidProperties'));
addpath(genpath('../../netLink/pipeFlow'));
addpath(genpath('../../netLink/conversionFactors'));



[dpVectorV,f_work,Z] = Beggs_Brill_dpds(s,pV,[],alpha,d,e,1,q_g_sc,q_o_sc,q_w_sc,rho_sc,s-1,s+1,T,T,hasSurfaceGas);

tic
[dpVectorV,f_work,Z] = Beggs_Brill_dpds(s,pV,[],alpha,d,e,1,q_g_sc,q_o_sc,q_w_sc,rho_sc,s-1,s+1,T,T,hasSurfaceGas);
vectorTime = toc


%[abs((dpVectorV-dpVector))*100./dpVectorV,dpVectorV,dpVector]

[v] = max(abs((dpVectorV-dpVector))*100./dpVectorV)



pScale  = 5*barsa;
qlScale =  5*meter^3/day;
qgScale = 100*(10*ft)^3/day;

pert = 1e-4;

qlP = qlScale*pert;
qgP = qgScale*pert;
pP =  pScale*pert;


gradScale = pScale./[pScale*ones(1,n),qgScale*ones(1,n),qlScale*ones(1,n),qlScale*ones(1,n)];


[dpVectorV] = Beggs_Brill_dpds(s,pV,[],alpha,d,e,1,q_g_sc,q_o_sc,q_w_sc,rho_sc,s-1,s+1,T,T,hasSurfaceGas);
tic
Jo =  (+Beggs_Brill_dpds(s,pV,[],alpha,d,e,1,q_g_sc,q_o_sc+qlP,q_w_sc,rho_sc,s-1,s+1,T,T,hasSurfaceGas)...
      -Beggs_Brill_dpds(s,pV,[],alpha,d,e,1,q_g_sc,q_o_sc-qlP,q_w_sc,rho_sc,s-1,s+1,T,T,hasSurfaceGas))/(2*qlP);

Jw =  (+Beggs_Brill_dpds(s,pV,[],alpha,d,e,1,q_g_sc,q_o_sc,q_w_sc+qlP,rho_sc,s-1,s+1,T,T,hasSurfaceGas)...
      -Beggs_Brill_dpds(s,pV,[],alpha,d,e,1,q_g_sc,q_o_sc,q_w_sc-qlP,rho_sc,s-1,s+1,T,T,hasSurfaceGas))/(2*qlP);

Jg =  (+Beggs_Brill_dpds(s,pV,[],alpha,d,e,1,q_g_sc+qgP,q_o_sc,q_w_sc,rho_sc,s-1,s+1,T,T,hasSurfaceGas)...
      -Beggs_Brill_dpds(s,pV,[],alpha,d,e,1,q_g_sc-qgP,q_o_sc,q_w_sc,rho_sc,s-1,s+1,T,T,hasSurfaceGas))/(2*qgP);
  
Jp =  (+Beggs_Brill_dpds(s,pV+pP,[],alpha,d,e,1,q_g_sc,q_o_sc,q_w_sc,rho_sc,s-1,s+1,T,T,hasSurfaceGas)...
      -Beggs_Brill_dpds(s,pV-pP,[],alpha,d,e,1,q_g_sc,q_o_sc,q_w_sc,rho_sc,s-1,s+1,T,T,hasSurfaceGas))/(2*pP);  
jacDFtime = toc


tic
[pVAD,q_g_scAD,q_o_scAD,q_w_scAD] =  initVariablesADI(pV,q_g_sc,q_o_sc,q_w_sc);
[dpVectorVAD2] = Beggs_Brill_dpds(s,pVAD,[],alpha,d,e,1,q_g_scAD,q_o_scAD,q_w_scAD,rho_sc,s-1,s+1,T,T,hasSurfaceGas);
jacADIT2 = toc
tic
[pVAD,q_g_scAD,q_o_scAD,q_w_scAD] =  initVariablesADI(pV,q_g_sc,q_o_sc,q_w_sc);
[dpVectorVAD] = Beggs_Brill_dpds(s,pVAD,[],alpha,d,e,1,q_g_scAD,q_o_scAD,q_w_scAD,rho_sc,s-1,s+1,T,T,hasSurfaceGas,f_work,Z);
jacADIT = toc

 
%adopts = admOptions('independents', [2,8,9,10],'dependents',1,'f', '-I /home/codas/Dropbox/PhD/mrst/remsoProject/netLink/pipeFlow/:/home/codas/Dropbox/PhD/mrst/remsoProject/netLink/fluidProperties/:/home/codas/Dropbox/PhD/mrst/remsoProject/netLink/conversionFactors');

%tic
%[JF] = admDiffVFor(@Beggs_Brill_dpds, 1, s,double(pVAD),[],alpha,d,e,1,double(q_g_scAD),double(q_o_scAD),double(q_w_scAD),rho_sc,s-1,s+1,T,T,hasSurfaceGas,f_work,Z, adopts);
%jacAdimatTime = toc


%max(max(abs(bsxfun(@rdivide,(cell2mat(dpVectorVAD.jac)-JF),gradScale))))
max(max(abs(bsxfun(@rdivide,full(cell2mat(dpVectorVAD.jac)-cell2mat(dpVectorVAD2.jac)),gradScale))))
max(max(abs(bsxfun(@rdivide,(cell2mat(dpVectorVAD.jac)-[diag(Jp),diag(Jg),diag(Jo),diag(Jw)]),gradScale))))


