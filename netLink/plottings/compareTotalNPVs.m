function [] = compareTotalNPVs(mFiles, obj, ss, times, scale)
%UNTITLED Compare NPV results from two different strategies    
    totalObj = cell(numel(mFiles),1);
    for i=1:numel(mFiles)          
        fileName = char(mFiles(i));       
        totalObj{i} = plotNPV(fileName, obj, ss, times, scale);
                 
        hold on;
        hold all;            
    end                  
    
    legendTxt = arrayfun(@(x1,x2) strcat('NPV: $', num2str(round(cell2mat(x2))), ', Strategy: ', cell2mat(x1)) ,mFiles, totalObj, 'UniformOutput', false);
    legend(strrep(legendTxt, '.mat', ''));
end


function [totalObj] = plotNPV(mFile, obj, ss, times, scale)
    barrelPrice = 30; %% in dollars
    load(mFile);
    totalPredictionSteps = getTotalPredictionSteps(ss);
    objs = zeros(1,totalPredictionSteps);
    for k = 1:totalPredictionSteps
        if k > 1
            [objs(k)] = callArroba(obj{k},{x{k-1},u{callArroba(ss.ci,{k})},v{k}}).*barrelPrice./scale;
        elseif k == 1
            [objs(k)] = callArroba(obj{k},{ss.state,u{callArroba(ss.ci,{k})},v{k}}).*barrelPrice./scale;
        else
            error('what?')
        end    
    end
    totalObj = -sum(objs);
    
    %plot results
    
    figure(1);   
    plot(times.steps(2:end), cumsum(-objs), '-x');
    %     plot(times.tPieceSteps, objAcum, '-x')
    ylabel('Accumulated NPV ($)');
    xlabel('time (day)');
    title('Cumulative Objective');   
end

