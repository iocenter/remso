function vert = newVertex(nV, sign, ws)
%%TODO: include name of vertices, use varargin e opt instead of parameters
%%list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% newVertex: creates a new connection point (vertex) in the graph  %
%%                                                                  %                
%% sign:  -1 for uphill flow, 1 for downhill flow                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
   vert = struct('id', nV, ...
                'name', '', ...
                'sign', sign, ...      % flow direction                 
                'Ein', [], ...         % set of edges entering node vertex
                'Eout', []);        % set of edges leaving node vertex
            
   if (nargin == 3) % wellSol is given   
       vert.name = ws(nV).name;       
   end
end
