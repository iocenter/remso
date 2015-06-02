function [adjRes,RHS] = solveFwdTransportSystem(G, S, W, rock, fluid, simRes, adjRes, obj, varargin)

% Find current time step (search for empty slots in adjRes)
% NOTE: actually curent time step +1

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



dt      = simRes(curStep).timeInterval * [-1 1]';

% Generate system matrix

numC    = G.cells.num;
PV      = G.cells.volumes.*rock.poro;
invDPV  = spdiags(1./PV, 0, numC, numC);

[mob, dmob] = mobilities(simRes(curStep).resSol, fluid);
Lt          = sum(mob, 2);

f           = bsxfun(@rdivide, mob, Lt);
Dfw_all     = (f(:,2).*dmob(:,1) - f(:,1).*dmob(:,2)) ./ Lt;   % Chain rule.
DDf         = spdiags(Dfw_all, 0, numC, numC);

%DDf     = spdiags( fluid.dfw(simRes(curStep).resSol), 0, numC, numC);


At      = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
    simRes(curStep).wellSol, 'Transpose', true);
systMat = speye(numC, numC) - dt * ( DDf * At * invDPV);   % system matrix
clear PV invDPV DDf At



ds =  - systMat'\ sRightSeed;
end

%--------------------------------------------------------------------------

function [mob, dmob] = mobilities(state, fluid)
%output/derivatives should be wrt s_w
mu = fluid.properties(state);
s  = fluid.saturation(state);
[kr{1:2}] = fluid.relperm(s, state);

%        \lambda_i in varargout{1}.
% (d/ds) \lambda_i in varargout{2}.  Returned iff requested.
%
mob = bsxfun(@rdivide, kr{1}, mu);
if nargout > 1
    dmob = bsxfun(@rdivide, kr{2}(:, [1 end]), mu);
    dmob(:, 2) = -dmob(:,2);
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
