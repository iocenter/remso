function vert = newVertex(nV, sign, tp, ws)
%%TODO: include name of vertices, use varargin e opt instead of parameters
%%list
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
                'name', '', ...
                'sign', sign, ...
                'type', tp, ...
                'pressure',0, ...
                'Ein', [], ...         % set of edges entering node vertex
                'Eout', [], ...        % set of edges leaving node vertex
                'qoV', 0, ...          % total of oil at a node vertex
                'qgV', 0, ...          % total of gas at a node vertex
                'qwV', 0, ...          % total of water at a node vertex      
                'flagStop', false);  % flag indicating the vertex in which the pressure back calculation should stop.
     
            
   if (nargin == 4) % wellSol is given   
       vert.name = ws(nV).name;
       vert.pressure = ws(nV).bhp/barsa;
       vert.qoV = ws(nV).qOs/(meter^3/day);
       vert.qgV = ws(nV).qGs/(meter^3/day);
       vert.qwV = ws(nV).qWs/(meter^3/day);
       
   end
end
