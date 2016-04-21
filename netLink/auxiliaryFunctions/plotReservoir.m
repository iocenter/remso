function [] = plotReservoir(reservoirP, W, varargin)
%%plotReservoir this function displays a saturation or pressure map for the
%%reservoir during the reservoir simulation.
    opt = struct('showPressures'       , false      , ...
                 'showSaturations'   , true, ...
                 'verbose', false);            

    opt = merge_options(opt, varargin{:});
    
    if opt.verbose
        mrstVerbose on;
    end

    G = reservoirP.G;
    initState = reservoirP.state;
    rock = reservoirP.rock;
    system = reservoirP.system;
    schedule = reservoirP.schedule;


    [sol, states, its] = runScheduleADI(initState, G, rock, system, schedule);

    clf;

    if opt.showPressures
        title(' Pressure Map');
    elseif opt.showSaturations
        title(' Saturation Map');
    end    

    colorbar;
    view(30,50);

    plotWell(G, W);
    view(30,50);

    for i=1:numel(states)
        if opt.showPressures
            plotCellData(G, states{i}.pressure./barsa) % pressure
        end
        
        if opt.showSaturations
            plotCellData(G, states{i}.s(:,2)) % saturation
        end
        pause(0.01);
    end
    
    
    
    [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
    clf;
    plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
%     plotCellData(G, reservoirP.state.pressure./barsa, j == round(G.cartDims(2)/2))
    plotCellData(G, reservoirP.state.pressure./barsa);
    % Plot the wells
    plotWell(G, W);
    view(30,50)
%     colorbar;
    
    [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
    clf;
    plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
    plotCellData(G, states{50}.pressure./barsa, j == round(G.cartDims(2)/2))
    % Plot the wells
    plotWell(G, W);
    view(30,50)
    
    
end