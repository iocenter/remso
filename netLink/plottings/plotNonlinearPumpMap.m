function [] = plotNonlinearPumpMap(dpFlows, dpPressures, times, qf_start, qf_end, numFlows, freq_start, freq_end, numFreq, nStages, fref, qmin_60, qmax_60)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plotting frequency bounds %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% TODO: receive netsol mock object and get densities from the network stream.
    oilDens = 897;
    waterDens = 1025.2;

    dp = zeros(numFlows, numFreq);
    flows = zeros(numFlows, numFreq);
    freq = zeros(numFreq, 1);
    interFlows = 1:numFlows;
    interFreq  = 1:numFreq;

    for j=interFreq
        f = freq_start+ (j-1).*(freq_end-freq_start)./(numFreq-1);
        freq(j)  =  f;
        for k=interFlows
            q = qf_start(j)+ (k-1).*(qf_end(j)-qf_start(j))./(numFlows-1);            
            flows(k,j) =  q./(meter^3/day);            
            if f == freq_end
                mixDens = 0.6*waterDens + 0.4*oilDens;
            else
                mixDens = oilDens;
            end
            
            dp(k,j) = calcDp(q, f, fref, nStages, 'mixDensity', mixDens);
            
        end
    end

    for k=interFreq
        plot(flows(:,k), dp(:,k), 'rx-');
        hold on;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plotting rate bounds %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numFreq = 10;
    interFreq  = 1:numFreq;
    qf_min = zeros(numFreq, 1);
    qf_max = zeros(numFreq, 1);
    incDens = 0.06;
    mixDens = oilDens;
    for j = interFreq(1:end)
        f1 = freq_start+ (j-1).*(freq_end-freq_start)./numFreq;
        f2 = freq_start+ (j).*(freq_end-freq_start)./numFreq;


        qf_min(j) = pump_rate(f1, qmin_60, 60);
        qf_min(j+1) = pump_rate(f2, qmin_60, 60);

        dp_min_1 = calcDp(qf_min(j),f1, fref, nStages, 'mixDensity',mixDens );
        dp_min_2 = calcDp(qf_min(j+1),f2, fref, nStages, 'mixDensity',mixDens);

        mixDens = (1-incDens*j)*(oilDens) + (incDens*j)*waterDens; %% increases water density up to 0.6 and decreases oil density down to 0.4

        line([qf_min(j)./(meter^3/day) qf_min(j+1)./(meter^3/day)], [dp_min_1,  dp_min_2], 'Marker','*','LineStyle','-', 'Color',[1 0 0]);
        hold on;

        qf_max(j) = pump_rate(f1, qmax_60, 60);
        qf_max(j+1) = pump_rate(f2, qmax_60, 60);

        dp_max_1 = calcDp(qf_max(j),f1, fref, nStages, 'mixDensity',mixDens );
        dp_max_2 = calcDp(qf_max(j+1),f2, fref, nStages, 'mixDensity',mixDens );


        line([qf_max(j)./(meter^3/day) qf_max(j+1)./(meter^3/day)], [dp_max_1,  dp_max_2], 'Marker','*','LineStyle','-', 'Color',[1 0 0]);
        hold on;

        qmedK1 = (qf_max(j)+qf_min(j))./2;
        qmedK2 = (qf_max(j+1)+qf_min(j+1))./2;

        dp_medK1 = calcDp(qmedK1,f1, fref, nStages, 'mixDensity',mixDens);
        dp_medK2 = calcDp(qmedK2,f2, fref, nStages, 'mixDensity',mixDens);

        line([qmedK1./(meter^3/day) qmedK2./(meter^3/day)], [dp_medK1,  dp_medK2],'LineStyle','--', 'Color',[0 0.64 0]);
    end

    hold on;             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plotting actual qf vs dp in the maps  %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    colorScale = [0 0 255]./255;
    for j=1:numel(times.steps(2:end))             
        colorScale = colorScale - [0 0 1]./255;
        plot(dpFlows(j)./(meter^3/day), -dpPressures(j)./barsa, 'Color',colorScale, 'Marker','*');
        hold on
    end
    
end