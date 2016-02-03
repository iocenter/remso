function [figN] = plotLinearPumpConstraints(flows, dp, times, netCst, figN, varargin)
    opt     = struct('extremePoints', []);
    opt     = merge_options(opt, varargin{:});
        
    if  ~isempty(opt.extremePoints)
        qminFmin = cell2mat(opt.extremePoints(1));
        qminFmax = cell2mat(opt.extremePoints(2));
        qmaxFmin = cell2mat(opt.extremePoints(3));
        qmaxFmax = cell2mat(opt.extremePoints(4));               
        
        for i=1:netCst
            figure(figN);
            figN = figN + 1;
            title(strcat('Pump Map: p ', int2str(i)));
            xlabel('Flow rate (sm3/day)');           
            ylabel('Pressure Difference (bar)');                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Plotting Linear Approx. of Pump Maps %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % [m1, n1] =  linearCoefficients(qminFmin, qmaxFmin);
            line([qminFmin(i,1) qmaxFmin(i,1)], [qminFmin(i,2),  qmaxFmin(i,2)], 'Color',[0 0 0]); % l1

            %  [m2, n2]=  linearCoefficients(qmaxFmin, qmaxFmax);
            line([qmaxFmin(i,1) qmaxFmax(i,1)], [qmaxFmin(i,2),  qmaxFmax(i,2)], 'Color',[0 0 0]); % l2

            %     [m3, n3] =  linearCoefficients(qminFmax, qmaxFmax);
            line([qminFmax(i,1) qmaxFmax(i,1)], [qminFmax(i,2),  qmaxFmax(i,2)], 'Color',[0 0 0]); % l3

            %    [m4, n4] =  linearCoefficients(qminFmin, qminFmax);
             line([qminFmin(i,1) qminFmax(i,1)], [qminFmin(i,2),  qminFmax(i,2)], 'Color',[0 0 0]); % l4                    
             
            hold all, hold on;  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Plotting Original (nonlinear) pump maps %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% TODO: plot nonlinear pump maps
            qf_start = [qminFmin(i,1); qminFmax(i,1)];
            qf_end = [qmaxFmin(i,1); qmaxFmax(i,1)];
            numFlows = 10;
            
            freq_start = 30;
            freq_end  = 90;
            numFreq  = 2;
            
            numStages = 80;
            fref = 60;
            
            qmin_60 = 15*(meter^3/day);
            qmax_60 = 200*(meter^3/day);    
            
            plotNonlinearPumpMap(qf_start, qf_end, numFlows, freq_start, freq_end, numFreq, numStages, fref, qmin_60, qmax_60);                        

            hold on;             
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Plotting actual qf vs dp in the maps  %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
                   
            for j=1:numel(times.steps(2:end))             
                plot(flows(i,j)./(meter^3/day), abs(dp(i,j))./barsa, 'bx-');
                hold on
            end
        end
        
    end
end


function [] = plotNonlinearPumpMap(qf_start, qf_end, numFlows, freq_start, freq_end, numFreq, nStages, fref, qmin_60, qmax_60)
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plotting frequency bounds %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        for i=interFlows
            q = qf_start(j)+ (i-1).*(qf_end(j)-qf_start(j))./(numFlows-1);

            flows(i,j) =  q;
            if f == freq_end
                mixDens = 0.6*waterDens + 0.4*oilDens;
            else
                mixDens = oilDens;
            end
                
            dp(i,j) = calcDp(q*(meter^3/day), f, fref, nStages, 'mixDensity', mixDens);

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
        
        mixDens = (1-incDens*j)*(oilDens) + (incDens*j)*waterDens; %% increases water density up to 0.6 and decreases oil density up to 0.4        
                       
        line([qf_min(j)./(meter^3/day) qf_min(j+1)./(meter^3/day)], [dp_min_1,  dp_min_2], 'Marker','*','LineStyle','-', 'Color',[1 0 0]);
        hold on;
        
        qf_max(j) = pump_rate(f1, qmax_60, 60);
        qf_max(j+1) = pump_rate(f2, qmax_60, 60);        

        dp_max_1 = calcDp(qf_max(j),f1, fref, nStages, 'mixDensity',mixDens );
        dp_max_2 = calcDp(qf_max(j+1),f2, fref, nStages, 'mixDensity',mixDens );
        
        
        line([qf_max(j)./(meter^3/day) qf_max(j+1)./(meter^3/day)], [dp_max_1,  dp_max_2], 'Marker','*','LineStyle','-', 'Color',[1 0 0]);
        hold on;
    end
end


















