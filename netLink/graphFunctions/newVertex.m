function vert = newVertex(nV, sign, tp, ws)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% newVertex: creates a new connection point (vertex) in the graph  %
%%                                                                  %                
%% sign:  -1 for uphill flow, 1 for downhill flow                   %
%%                                                                  %
%% type:  0 for connection vertex,                                  %                
%%       -1 for production ending vertex,                           %
%%       -2 for production starting,                                % 
%%        1 for injection starting vertex,                          % 
%%        2 for injection ending vertex                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
   vert = struct('id', nV, ...
                'sign', sign, ...
                'type', tp, ...
                'pressure',0, ...
                'Ein', [], ...          % set of edges entering node vertex
                'Eout', [], ...         % set of edges leaving node vertex
                'qoV', 0, ...          % total of oil at a node vertex
                'qgV', 0, ...          % total of gas at a node vertex
                'qwV', 0);            % total of water at a node vertx      
     
            
   if (nargin == 4) % wellSol is given   
       vert.pressure = ws(nV).bhp*1e-05;
       vert.qoV = ws(nV).qOs*day;
       vert.qgV = ws(nV).qGs*day;
       vert.qwV = ws(nV).qWs*day;
%        vert = struct('id', nV,'sign', sign,'type', tp,'pressure', ws(nV).bhp, 'qoV', ws(nV).qOs, 'qwV', ws(nV).qWs, 'qgV', ws(nV).qGs);               
   end
end
