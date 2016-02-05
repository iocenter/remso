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
end