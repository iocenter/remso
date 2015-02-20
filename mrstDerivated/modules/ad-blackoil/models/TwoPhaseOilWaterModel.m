classdef TwoPhaseOilWaterModel < ReservoirModel
    % Two phase oil/water system without dissolution

%{
Changes by Codas

model.scaling
function varargout = toMRSTStates(model,stateVector)          
function varargout = toStateVector(model,state)

%}
    properties

    end
    
    methods
        function model = TwoPhaseOilWaterModel(G, rock, fluid, varargin)
            
            model = model@ReservoirModel(G, rock, fluid);
            
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
            
            model = merge_options(model, varargin{:});
            
            % Setup operators
            model = model.setupOperators(G, rock, 'deck', model.inputdata);
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWater(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        
        function varargout = toMRSTStates(model,stateVector)
            
            nx = model.G.cells.num;
            
            stateMrst.pressure = stateVector(1:nx)*model.scaling.p;
            stateMrst.s = [stateVector(nx+1:end),1-stateVector(nx+1:end)]*model.scaling.s;
            
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
        
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            % Parent class handles almost everything for us
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
            
            % Update wells based on black oil specific properties
            saturations = model.saturationVarNames;
            wi = strcmpi(saturations, 'sw');
            oi = strcmpi(saturations, 'so');
            gi = strcmpi(saturations, 'sg');

            W = drivingForces.Wells;
            state.wellSol = assignWellValuesFromControl(model, state.wellSol, W, wi, oi, gi);

        end
    end
end
