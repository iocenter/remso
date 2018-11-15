function netSol = prodNetwork(wellSol, varargin)
%%
%  Creates an object containing the topology of a production/injection 
%  network that is compatible the wellSol object.
%
% SYNOPSIS:
%  [u,x,v,f,xd,M,simVars] = prodNetwork(wellSol, ...)
% PARAMETERS:
%   wellSol - Algebraic variables for the wells containing rates and
%   pressures at a give time interval%
%
% RETURNS:
%
%   netSol - Network object containing the following attributes:
%               V: All vertices of the network
%            Vsrc: Set of source vertices
%            Vsnk: Set of sink vertices
%              Vw: Set of well vertices (both injectors and producers)
%          VwProd: Set of production well vertices
%           VwInj: Set of injection well vertices
%               E: Set of edges
%            Eeqp: Set of 'special' edges, i.e., edges denoting equipment
%                  in the network such as ESPs or Chokes
%               M: Flow Matrix, which is used to compute the flows in the
%                  edges from the well flow rates with a linear operation
%              qo: Oil flow rates in the edges
%              qw: Water flow rates in the edges
%              qg: Gas flow rates in the edges
%              pV: Pressures in the vertices
%              boundaryCond: boudary conditios, e.g., constant pressure in
%                            the inlet of the separator at the surface 
%                            facilities.
%
%{

Copyright 2015-2018, Thiago Lima Silva

REMSO is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

REMSO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with REMSO.  If not, see <http://www.gnu.org/licenses/>.

%} 
    opt = struct('simpleNetwork',false, 'toyNetwork',false, 'eirikNetwork', false, 'espNetwork', false, 'satelliteWellsNetwork', false, 'withPumps', false);
    opt = merge_options(opt, varargin{:});

    netSol = initNetSolLocal(wellSol);       
    
    if opt.simpleNetwork
        netSol = createSimpleNetwork(netSol);  
    elseif opt.toyNetwork    
        netSol = createToyNetwork(netSol);
    elseif opt.eirikNetwork
        netSol = createEirikNetwork(netSol);
    elseif opt.espNetwork
        netSol = createESPNetwork(netSol, 'withPumps', opt.withPumps);
    elseif opt.satelliteWellsNetwork
        netSol = createSatelliteWellsNetwork(netSol);
    end
    
    % flow matrix
    netSol.M = createFlowMatrix(netSol);
    
    % init flow vectors
    netSol.qo = zeros(numel(netSol.E),1);
    netSol.qw = zeros(numel(netSol.E),1);
    netSol.qg = zeros(numel(netSol.E),1);        
    
    % init pressure vector
    netSol.pV = zeros(numel(netSol.V),1);
    netSol.pV(setdiff(netSol.Vsnk, netSol.VwInj)) = netSol.boundaryCond; % network boundary condition
end  

