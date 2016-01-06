function [] = plotBO(fluid )

if isempty(fluid)
    deck = readEclipseDeck('odeh_adi.data');
    
    deck = convertDeckUnits(deck);
    
    fluid = initDeckADIFluid(deck);
end


fVars = functions(fluid.rsSat);
prange = [min(fVars.workspace{1}.pvto{1}.data(:,1)),max(fVars.workspace{1}.pvto{1}.data(:,1))];

ps = (prange(1):(prange(2)-prange(1))/1000:prange(2))';

rsSat = fluid.rsSat(ps);


psSAT = (prange(1):(prange(2)-prange(1))/10:prange(2)-(prange(2)-prange(1))/10)';

rsSatV = fluid.rsSat(psSAT);

figure(1) ; clf; hold all;


flag = true(size(rsSat));
bOS = fluid.bO(ps,rsSat,flag);
bOU = fluid.bO(ps,rsSat,~flag);

plot(ps,bOS,'-.',ps,bOU,':')

for k = 1:numel(rsSatV)
    
    psu = (psSAT(k):(prange(2)-prange(1))/10000:prange(2))';
    rs = ones(size(psu))*rsSatV(k);
    flag = false(size(psu));
    
    bO = fluid.bO(psu,rs,flag);
    
    plot(psu,bO)
    
end


title('bO')
legend('fluid.bO(p,rsSat(p),true)','fluid.bO(p,rsSat(p),false)','Others are Undersatured for fixed rs')
