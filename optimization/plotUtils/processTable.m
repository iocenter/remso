function [it,obj,lag,du,eq,ineq,rho,tau,eff,qpi,norm1QpViol,cond ] = processTable(table)



% remove headers
k = 1;
while k <= size(table,1)
    while strcmp(table{k,1}(end),' ') && length(table{k,1}) >= 2
        table{k,1} = table{k,1}(1:end-1);
    end
    if strcmp(table(k,1),'it')
        table = [table(1:k-1,:);table(k+1:end,:)];
    elseif strcmp(table{k,1}(end),'*')
        table{k,1} = table{k,1}(1:end-1);
    end

  
    k = k+1;
end

%  getObjectivePerItartion

itPos = cellfun(@(x)~isempty(toNum(x)),table(:,1));
itFPos = [itPos(2:end);true];
it = cellfun(@(x)toNum(x),table(itPos,1));

obj = cat(1,table{itPos,2},table{end,2});


lag = cellfun(@(x)toNum(x),table(itPos,3));
du = cellfun(@(x)x,table(itFPos,4));
eq = cat(1,table{itPos,5},table{end,5});
ineq = cat(1,table{itPos,6},table{end,6});

rho = cellfun(@(x)toNum(x),table(itPos,7));
tau = cellfun(@(x)toNum(x),table(itPos,8));
eff = cat(1,table{1,12},table{~itPos,12});
qpi  = cellfun(@(x)toNum(x),table(itPos,13));
norm1QpViol = cellfun(@(x)toNum(x),table(itPos,14)); 
cond = cellfun(@(x)toNum(x),table(itPos,15)); 




end

function n = toNum(x)

    if isnumeric(x)
        n = x;
    else
        n = str2num(x);
    end

end