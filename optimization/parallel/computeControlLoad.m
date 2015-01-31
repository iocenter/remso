function [ uLoad,uStart ] = computeControlLoad(ci,totalControlSteps,totalPredictionSteps)


uStart = inf(1,totalControlSteps);

for k = 1:totalPredictionSteps
    uk = callArroba(ci,{k});
    uStart(uk) = min(uStart(uk),k);

end

uLoad = zeros(1,totalControlSteps);

for k = 1:length(uLoad)
    uLoad(k) = (totalPredictionSteps-uStart(k)+1);
end


end

    