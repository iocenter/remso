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
   if (nargin == 4)
       vert = struct('id', nV,'sign', sign,'type', tp,'pressure', ws(nV).bhp, 'qoV', ws(nV).qOs, 'qwV', ws(nV).qWs, 'qgV', ws(nV).qGs);               
   else
       vert = struct('id', nV, 'sign', sign, 'type', tp, 'pressure', 0, 'qoV', 0, 'qwV', 0, 'qgV', 0);        
   end       
end
