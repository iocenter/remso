function objADI = simpleNPVADI(G, S, W, rock, fluid, simRes, ...
                         schedule, controls, varargin)
%Simple net-present-value function - no discount factor
%
% SYNOPSIS:
%   obj = simpleNPV(G, S, W, rock, fluid, simRes, schedule, controls)
%
% DESCRIPTION:
%   Computes value of objective function for given simulation, and partial
%   derivatives of variables if varargin > 6
%
% PARAMETERS:
%   simRes      -
%
% RETURNS:
%   obj         - structure with fields
%
% SEE ALSO:

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
%{
inspired in simpleNPV, using ADI and enabling leftSeed
adding scaling and sign
%}


opt     = struct('OilPrice'              , 100   , ...
                 'WaterProductionCost'   ,  10 , ...
                 'WaterInjectionCost'    ,  10 , ...
                 'RelativeDiscountFactor', 0  , ...
                 'ComputePartials',         [],...
                 'leftSeed',[],...
                 'sign',1,...
                 'scale',1);
opt     = merge_options(opt, varargin{:});

ro      = opt.OilPrice            / stb;
rw      = opt.WaterProductionCost / stb;
ri      = opt.WaterInjectionCost  / stb;
d       = opt.RelativeDiscountFactor;

%-----------------------------------------------
computePartials = opt.ComputePartials;
if isempty(computePartials)
    computePartials  = (nargin > 6);
end
numSteps = numel(simRes);
val      = 0;
dval = zeros(1,numSteps);
objADI = cell(1,numSteps);
partials = repmat( struct('v', [], 'p', [], 'pi', [], ...
                          's', [], 'u', []), [numSteps 1] );

totTime  = max( [simRes.timeInterval] );

mob = cell([1, 1 + 2*double(computePartials)]);

for step = 1 : numSteps,
    resSol  = simRes(step).resSol;
    wellSol = simRes(step).wellSol;
    int     = simRes(step).timeInterval;
    dt      = int(2) - int(1);
    dFac    = (1+d)^(-int(2)/totTime);
    dFac = dFac * opt.sign * opt.scale; % lets just add the scaling here

    [wellRates, rateSigns] = getRates(W, wellSol);
    wellCells = vertcat( W.cells );
    wellSats  = resSol.s( wellCells );

    [mob{:}] = mobilities(struct('s', wellSats), fluid);

    Lt  = sum(mob{1}, 2);
    f   = bsxfun(@rdivide, mob{1}, Lt);

    f_w = f(:,1);
    f_o = 1 - f_w;

    injInx  = (rateSigns > 0);
    prodInx = (rateSigns < 0);

    % Objective value:

    dval(step)   =  dt*dFac*( - sum(  wellRates(injInx)                )*ri ...
                            - sum( -wellRates(prodInx).*f_w(prodInx) )*rw ...
                            + sum( -wellRates(prodInx).*f_o(prodInx) )*ro );
    val = val+dval(step);
    if computePartials,
        numCF    = size(G.cells.faces, 1);
        numC     = G.cells.num;
        numF     = G.faces.num;
        numU     = numel(controls.well);

        partials(step).v   = zeros(1, numCF);
        partials(step).p   = zeros(1, numC);
        partials(step).pi  = zeros(1, numF);

        partials(step).q_w =  - dt*dFac*ri*injInx' ...
                              + dt*dFac*rw*( prodInx.*f_w )' ...
                              - dt*dFac*ro*( prodInx.*f_o )';

        df  = (f(:,2).*mob{2}(:,1) - f(:,1).*mob{2}(:,2)) ./ Lt;   % Chain rule.
        d2f = (f(:,2).*mob{3}(:,1) - f(:,1).*mob{3}(:,2) - ...
               2 .* df .* sum(mob{2}, 2)) ./ Lt;

        Df_w = df(:,1);
        Df_o = - Df_w;

        ds   = zeros(1, numC);
        ds( wellCells(prodInx) )  =  dt*dFac*rw*( wellRates(prodInx).*Df_w(prodInx) )' ...
                                    -dt*dFac*ro*( wellRates(prodInx).*Df_o(prodInx) )';

        partials(step).s = ds;
        partials(step).u = zeros(1, numU);


        objADI{step} = ADI(dval,{partials(step).s,partials(step).q_w});
        
        if size(opt.leftSeed,2) > 0;
            objADI{step}.jac = cellfun(@(x)opt.leftSeed*x,objADI{step}.jac,'UniformOutput',false);
        end

        
    else
        objADI{step} = dval;
    end
end


end

%--------------------------------------------------------------------------

function [wellRates, rateSigns] = getRates(W, wellSol)
   wellRates = vertcat(wellSol.flux);
   wellSys   = [ W.S ];
   Dw        = blkdiag( wellSys.D );
   wellSigns = ones( numel(W), 1 );
   totRates  = Dw'*wellRates;
   wellSigns( totRates < 0 ) = -1;

   for k = 1:numel(W),
      if ~isempty(W(k).sign) % override if sign is given expl
         wellSigns(k) = W(k).sign;
      end
   end

   rateSigns = Dw*wellSigns;
end

%--------------------------------------------------------------------------

function [mob, dmob, dmob2] = mobilities(state, fluid)
   %output/derivatives should be wrt s_w
   mu = fluid.properties(state);
   s  = fluid.saturation(state);
   [kr{1:nargout}] = fluid.relperm(s, state);

   %        \lambda_i in varargout{1}.
   % (d/ds) \lambda_i in varargout{2}.  Returned iff requested.
   %
   mob = bsxfun(@rdivide, kr{1}, mu);
   if nargout > 1
       dmob = bsxfun(@rdivide, kr{2}(:, [1 end]), mu);
       dmob(:, 2) = -dmob(:,2);
   end
    if nargout > 2
       dmob2 = bsxfun(@rdivide, kr{3}(:, [1 end]), mu);
       %dmob2(:, 2) = -dmob(:,2);
   end
   %kr = cellfun(@(x)x(:,[1 end]), kr, 'UniformOutput', false);
   %varargout = cellfun(@(n) bsxfun(@rdivide, n, mu), kr, ...
   %                    'UniformOutput', false);
end
%--------------------------------------------------------------------------

% function varargout = mobilities(state, fluid)
%    mu = fluid.properties(state);
%    s  = fluid.saturation(state);
%    [kr{1:nargout}] = fluid.relperm(s, state);
%
%    %        \lambda_i in varargout{1}.
%    % (d/ds) \lambda_i in varargout{2}.  Returned iff requested.
%    %
%    kr = cellfun(@(x)x(:,[1 end]), kr, 'UniformOutput', false);
%    varargout = cellfun(@(n) bsxfun(@rdivide, n, mu), kr, ...
%                        'UniformOutput', false);
% end

%{
function [wellRates, rateSigns] = getRates(W, wellSol)
wellRates   = vertcat(wellSol.flux);
wellSys     = [ W.S ];
Dw    = blkdiag( wellSys.D );
if isfield(W, 'sign')
    wellSigns = vertcat(W.sign);
else
    wellSigns    = ones( numel(W), 1 );
    totRates = Dw'*wellRates;
    wellSigns( totRates < 0 ) = -1;
end
%}
