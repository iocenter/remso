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
            title(strcat('Linear Pump Map: Well', int2str(i)));
            xlabel('flow (sm3/day)');           
            ylabel('Dp (bar)');                        
            
            % [m1, n1] =  linearCoefficients(qminFmin, qmaxFmin);
            line([qminFmin(i,1) qmaxFmin(i,1)], [qminFmin(i,2),  qmaxFmin(i,2)]); % l1

            %  [m2, n2]=  linearCoefficients(qmaxFmin, qmaxFmax);
            line([qmaxFmin(i,1) qmaxFmax(i,1)], [qmaxFmin(i,2),  qmaxFmax(i,2)]); % l2

            %     [m3, n3] =  linearCoefficients(qminFmax, qmaxFmax);
            line([qminFmax(i,1) qmaxFmax(i,1)], [qminFmax(i,2),  qmaxFmax(i,2)]); % l3

            %    [m4, n4] =  linearCoefficients(qminFmin, qminFmax);
             line([qminFmin(i,1) qminFmax(i,1)], [qminFmin(i,2),  qminFmax(i,2)]); % l4                    
            hold on, hold all;
            for j=1:numel(times.steps(2:end))             
                plot(flows(i,j)./(meter^3/day), abs(dp(i,j))./barsa, 'rx-');
                hold on
            end
        end
        
    end
end