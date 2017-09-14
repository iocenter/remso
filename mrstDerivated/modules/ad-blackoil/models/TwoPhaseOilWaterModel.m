classdef TwoPhaseOilWaterModel < ThreePhaseBlackOilModel
% Two phase oil/water system without dissolution
%{
Changes from https://github.com/iocenter/remso/ 4f8fa54a92c14b117445ad54e9e5dd3a0e47a7f5

model.scaling
function varargout = toMRSTStates(model,stateVector)          
function varargout = toStateVector(model,state)

%}
properties

end

methods
    function model = TwoPhaseOilWaterModel(G, rock, fluid, varargin)
        model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});

        % This is the model parameters for oil/water
        model.oil = true;
        model.gas = false;
        model.water = true;

        % Blackoil -> use CNV style convergence 
        model.useCNVConvergence = true;

        model.saturationVarNames = {'sw', 'so'};
        model.wellVarNames = {'qWs', 'qOs', 'bhp'};

        model.scaling = struct ('p',5*barsa,...
            's',0.01,...
            'qWs',10*meter^3/day,...
            'qOs',10*meter^3/day,...
            'bhp',5*barsa);
            
        % these options are processed on the ReservoirModel
        % included here to prevent a warning, opt will not be used
        opt = struct('allowWellSignChange', false, ...
                     'allowCrossflow', true, ...
                     'allowControlSwitching', true);
                     
        [opt, varargin] = merge_options(opt, varargin{:});             
            
        model = merge_options(model, varargin{:});
    end
    
    % --------------------------------------------------------------------%
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        [problem, state] = equationsOilWater(state0, state, model,...
                        dt, ...
                        drivingForces,...
                        varargin{:});

    end
    function varargout = toMRSTStates(model,stateVector)
        
        nx = model.G.cells.num;
        
        stateMrst.pressure = stateVector(1:nx)*model.scaling.p;
        sW = stateVector(nx+1:end)*model.scaling.s;
        stateMrst.s = [sW,1-sW];
        
        varargout{1} = stateMrst;
        
        if nargout >=2
            varargout{2} = sparse(1:2*nx,1:2*nx,[....
                ones(nx,1)*model.scaling.p;...
                ones(nx,1)*model.scaling.s]);
end
end

    function [varargout] = toStateVector(model,state)
        
        nx = model.G.cells.num;
        
        varargout{1} = [state.pressure/model.scaling.p;
            state.s(:,1)/model.scaling.s];
        
        if nargout >=2
            varargout{2} = sparse(1:2*nx,1:2*nx,[...
                ones(nx,1)/model.scaling.p;...
                ones(nx,1)/model.scaling.s]);
        end
        
    end
end
end
%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
