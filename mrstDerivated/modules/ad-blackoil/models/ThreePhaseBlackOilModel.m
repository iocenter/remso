classdef ThreePhaseBlackOilModel < ReservoirModel
% Three phase with optional dissolved gas and vaporized oil
%{
Changes from https://github.com/iocenter/remso/ 4f8fa54a92c14b117445ad54e9e5dd3a0e47a7f5

model.scaling
function varargout = toMRSTStates(model,stateVector)
function varargout = toStateVector(model,state)
testGetEquation and related functions for debuging --> Requires Adimat

%}
properties
    % Flag deciding if gas can be dissolved into the oil phase
    disgas
    % Flag deciding if oil can be vaporized into the gas phase
    vapoil

    % Maximum relative Rs/Rv increment
    drsMaxRel
    % Maximum absolute Rs/Rv increment
    drsMaxAbs
end

methods
    function model = ThreePhaseBlackOilModel(G, rock, fluid, varargin)
        model = model@ReservoirModel(G, rock, fluid, varargin{:});

        % Typical black oil is disgas / dead oil, but all combinations
        % are supported
        model.vapoil = false;
        model.disgas = false;

        % Max increments
        model.drsMaxAbs = inf;
        model.drsMaxRel = inf;

        % Blackoil -> use CNV style convergence 
        model.useCNVConvergence = true;

        % All phases are present
        model.oil = true;
        model.gas = true;
        model.water = true;
        model.saturationVarNames = {'sw', 'so', 'sg'};

        model.scaling = struct ('p',5*barsa,...
            's',0.01,...
            'rGH',0.01,...
            'qWs',10*meter^3/day,...
            'qOs',10*meter^3/day,...
            'qGs',1000*meter^3/day,...
            'bhp',5*barsa);
            
        % these options are processed on the ReservoirModel
        % included here to prevent a warning, opt will not be used
        opt = struct('allowWellSignChange', false, ...
                     'allowCrossflow', true, ...
                     'allowControlSwitching', true);
                     
        [opt, varargin] = merge_options(opt, varargin{:});               
            
        model = merge_options(model, varargin{:});

        d = model.inputdata;
        if ~isempty(d)
            % Assume ECL-style input deck, as this is the only
            % supported format at the moment.
            if isfield(d, 'RUNSPEC')
                if isfield(d.RUNSPEC, 'VAPOIL')
                    model.vapoil = d.RUNSPEC.VAPOIL;
                end
                if isfield(d.RUNSPEC, 'DISGAS')
                    model.disgas = d.RUNSPEC.DISGAS;
                end
            else
                error('Unknown dataset format!')
            end
        end
    end
    
    % --------------------------------------------------------------------%
    function [fn, index] = getVariableField(model, name)
        switch(lower(name))
            case {'rs', 'rv'}
                % RS and RV for gas dissolving into the oil phase and oil
                % components vaporizing into the gas phase respectively.
                fn = lower(name);
                index = 1;
            otherwise
                % Basic phases are known to the base class
                [fn, index] = getVariableField@ReservoirModel(model, name);
        end
    end
    
    % --------------------------------------------------------------------%
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        [problem, state] = equationsBlackOil(state0, state, model, dt, ...
                        drivingForces, varargin{:});

    end

    % --------------------------------------------------------------------%
    function state = validateState(model, state)
        % Check parent class
        state = validateState@ReservoirModel(model, state);
        nc = model.G.cells.num;
        if model.disgas
            % RS must be supplied for all cells. This may cause an error.
            model.checkProperty(state, 'rs', nc, 1);
        else
            % RS does not really matter. Assign single value.
            fn = model.getVariableField('rs');
            if ~isfield(state, fn)
                dispif(model.verbose, ...
                    ['Missing field "', fn, '" added since disgas is not enabled.\n']);
                state.(fn) = 0;
            end
            clear fn
        end
        if model.vapoil
            % RV must be supplied for all cells. This may cause an error.
            model.checkProperty(state, 'rv', nc, 1);
        else
            % RS does not really matter. Assign single value.
            fn = model.getVariableField('rv');
            if ~isfield(state, fn)
                dispif(model.verbose, ...
                    ['Missing field "', fn, '" added since vapoil is not enabled.\n']);
                state.(fn) = 0;
            end
            clear fn
        end
    end

    % --------------------------------------------------------------------%
    function [state, report] = updateState(model, state, problem, dx, drivingForces)
        saturations = lower(model.saturationVarNames);
        wi = strcmpi(saturations, 'sw');
        oi = strcmpi(saturations, 'so');
        gi = strcmpi(saturations, 'sg');

        vars = problem.primaryVariables;
        removed = false(size(vars));
        if model.disgas || model.vapoil
            % The VO model is a bit complicated, handle this part
            % explicitly.
            state0 = state;

            state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
            state = model.capProperty(state, 'pressure', model.minimumPressure, model.maximumPressure);

            [vars, ix] = model.stripVars(vars, 'pressure');
            removed(~removed) = removed(~removed) | ix;

            % Black oil with dissolution
            so = model.getProp(state, 'so');
            sw = model.getProp(state, 'sw');
            sg = model.getProp(state, 'sg');

            % Magic status flag, see inside for doc
            st = model.getCellStatusVO(state0, so, sw, sg);

            dr = model.getIncrement(dx, problem, 'x');
            dsw = model.getIncrement(dx, problem, 'sw');
            % Interpretation of "gas" phase varies from cell to cell, remove
            % everything that isn't sG updates
            dsg = st{3}.*dr - st{2}.*dsw;

            if model.disgas
                state = model.updateStateFromIncrement(state, st{1}.*dr, problem, ...
                                                       'rs', model.drsMaxRel, model.drsMaxAbs);
            end

            if model.vapoil
                state = model.updateStateFromIncrement(state, st{2}.*dr, problem, ...
                                                       'rv', model.drsMaxRel, model.drsMaxAbs);
            end

            dso = -(dsg + dsw);

            ds = zeros(numel(so), numel(saturations));
            ds(:, wi) = dsw;
            ds(:, oi) = dso;
            ds(:, gi) = dsg;

            state = model.updateStateFromIncrement(state, ds, problem, 's', model.dsMaxRel, model.dsMaxAbs);
            % We should *NOT* be solving for oil saturation for this to make sense
            assert(~any(strcmpi(vars, 'so')));
            state = computeFlashBlackOil(state, state0, model, st);
            state.s  = bsxfun(@rdivide, state.s, sum(state.s, 2));

            %  We have explicitly dealt with rs/rv properties, remove from list
            %  meant for autoupdate.
            [vars, ix] = model.stripVars(vars, {'sw', 'so', 'sg', 'rs', 'rv', 'x'});
            removed(~removed) = removed(~removed) | ix;

        end

        % We may have solved for a bunch of variables already if we had
        % disgas / vapoil enabled, so we remove these from the
        % increment and the linearized problem before passing them onto
        % the generic reservoir update function.
        problem.primaryVariables = vars;
        dx(removed) = [];

        % Parent class handles almost everything for us
        [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);

        % Handle the directly assigned values (i.e. can be deduced directly from
        % the well controls. This is black oil specific.
        W = drivingForces.W;
        state.wellSol = assignWellValuesFromControl(model, state.wellSol, W, wi, oi, gi);
    end
    
    % --------------------------------------------------------------------%
    function scaling = getScalingFactorsCPR(model, problem, names)
        % Get approximate, impes-like pressure scaling factors
        nNames = numel(names);
        
        scaling = cell(nNames, 1);
        handled = false(nNames, 1);
        
        % Take averaged pressure for scaling factors
        state = problem.state;
        fluid = model.fluid;
        p = mean(state.pressure);
        
        for iter = 1:nNames
            name = lower(names{iter});
            switch name
                case 'oil'
                    if model.disgas
                       rs = fluid.rsSat(p);
                       bO = fluid.bO(p, rs, true);
                    else
                       bO = fluid.bO(p);
                    end
                    s = 1./bO;
                case 'water'
                    bW = fluid.bW(p);
                    s = 1./bW;
                case 'gas'
                    if model.vapoil
                        rv = fluid.rvSat(p);
                        bG = fluid.bG(p, rv, true);
                    elseif model.gas
                        bG = fluid.bG(p);
                    end
                    s = 1./bG;
                otherwise
                    continue
            end
            sub = strcmpi(problem.equationNames, name);
            
            scaling{iter} = s;
            handled(sub) = true;
        end
        if ~all(handled)
            % Get rest of scaling factors from parent class
            other = getScalingFactorsCPR@ReservoirModel(model, problem, names(~handled));
            [scaling{~handled}] = other{:};
        end
    end
    function varargout = toMRSTStates(model,stateVector)
    
        partials = nargout >=2 ;
        
        nx = model.G.cells.num;
        
        p = stateVector(1:nx)*model.scaling.p;
        sW = stateVector(nx+1:2*nx)*model.scaling.s;
        rGH = stateVector(2*nx+1:end)*model.scaling.rGH;
        
        [ stateMrst,Jac ] = statePsWrGH2stateMRST( p,sW,rGH,model.fluid,model.disgas,model.vapoil,'partials',partials);
        
        varargout{1} = stateMrst;
        
        if partials
            Jac{1} = Jac{1}*model.scaling.p;
            Jac{2} = Jac{2}*model.scaling.s;
            Jac{3} = Jac{3}*model.scaling.rGH;
            
            varargout{2} = cat(2,Jac{:});
        end
    end
    function [varargout] = toStateVector(model,state)
        
        partials = nargout >=2 ;
        
        [ p,sW,rGH ] = stateMrst2statePsWrGH(state,model.fluid,model.disgas,model.vapoil,'partials',partials);
        
        stateVector = [p  /model.scaling.p;
            sW /model.scaling.s;
            rGH/model.scaling.rGH];
        
        
        varargout{1} = double(stateVector);
        
        if partials
            
            stateVector = cat(stateVector);
            varargout{2} = stateVector.jac{1};
            
        end
    end    
    function st = getCellStatusVO(model, state, sO, sW, sG)
        status = [];
        if isfield(state, 'status')
            status = state.status;
        end
        st = getCellStatusVO(sO, sW, sG, 'status', status, 'vapoil', ...
                                 model.vapoil, 'disgas', model.disgas);
    end
    
    function [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, ...
                                                          status, sO, x, rs, ...
                                                          rv, pressure)
        [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model.fluid, ...
                                                          status, sO, x, rs, ...
                                                          rv, pressure, model.disgas, model.vapoil);
    end
    
    function [e] = testGetEquation(model, state0, state, dt, drivingForces,varargin)
        
        opt.fdStep = sqrt(eps)*100;
        
        [opt,otherOpts] = merge_options(opt,varargin{:});
        
        e = 0;
        
        
        [residualScale] = getResidualScale(model,state,dt);
        
        
        % Test Jacobian w.r.t state0
        state0Vector = model.toStateVector(state0);
        
        
        [V,Jac0] = testGetEquationState0(model, state0Vector, state, dt, drivingForces,residualScale,varargin{:});
        f0 = @(x) testGetEquationState0(model, x, state, dt, drivingForces,residualScale,varargin{:});
        
        ADopts = admOptions();
        ADopts.JPattern = Jac0;
        ADopts.fdStep = opt.fdStep;
        Jac02 = admDiffFD(f0, 1, state0Vector, ADopts);
        Jac02 = sparse(Jac02);
        
        
        e = [e max(max(abs(Jac0-Jac02)))];
        
        % Test Jacobian w.r.t state
        stateVector = model.toStateVector(state);
        [V,Jac] = testGetEquationState(model, state0, stateVector,state.wellSol, dt, drivingForces,residualScale,varargin{:});
        f0 = @(x) testGetEquationState(model, state0, x,state.wellSol, dt, drivingForces,residualScale,varargin{:});
        
        ADopts = admOptions();
        ADopts.JPattern = Jac;
        %size(admColorSeed(Jac, ADopts),2); % this number is proportional to the callings to be done
        ADopts.fdStep = opt.fdStep;
        Jac2 = admDiffFD(f0, 1, stateVector, ADopts);
        Jac2 = sparse(Jac2);
        
        
        e = [e max(max(abs(Jac-Jac2)))];
        
        
        % Test Jacobian w.r.t wellSol
        wellSolVector = model.toWellSolVector(state.wellSol);
        [Vw,JacW] = testGetEquationWell(model, state0, state,wellSolVector, dt, drivingForces,residualScale,varargin{:});
        f0 = @(x) testGetEquationWell(model, state0, state,x, dt, drivingForces,residualScale,varargin{:});
        
        ADopts = admOptions();
        ADopts.JPattern = JacW;
        ADopts.fdStep = opt.fdStep;
        JacW2 = admDiffFD(f0, 1, wellSolVector, ADopts);
        JacW2 = sparse(JacW2);
        
        
        
        
        e = [e max(max(abs(JacW-JacW2)))];
    
end
    function [varargout] = testGetEquationState0(model, state0Vector, state, dt, drivingForces,residualScale,varargin)
        
        
        if nargout == 1
            varargin = [varargin,{'resOnly',true}];
            
            [state0Mrst] = model.toMRSTStates(state0Vector);
            [problem] = getEquations(model, state0Mrst, state, dt, drivingForces,varargin{:},'reverseMode',true);
            
            problem.equations{1} = problem.equations{1}./residualScale.o;
            problem.equations{2} = problem.equations{2}./residualScale.w;
            problem.equations{3} = problem.equations{3}./residualScale.g;
            problem.equations{4} = problem.equations{4}/model.scaling.qWs;
            problem.equations{5} = problem.equations{5}/model.scaling.qOs;
            problem.equations{6} = problem.equations{6}/model.scaling.qGs;
            problem.equations{7} = problem.equations{7}/model.scaling.bhp;
            
            
            
            Val = vertcat(problem.equations{:});
            
            varargout{1} = Val;
            
        else
            varargin = [varargin,{'resOnly',false}];
            
            [state0Mrst,JacX0] = model.toMRSTStates(state0Vector);
            [problem] = getEquations(model, state0Mrst, state, dt, drivingForces,varargin{:},'reverseMode',true);
            
            
            problem.equations{1} = problem.equations{1}./residualScale.o;
            problem.equations{2} = problem.equations{2}./residualScale.w;
            problem.equations{3} = problem.equations{3}./residualScale.g;
            problem.equations{4} = problem.equations{4}/model.scaling.qWs;
            problem.equations{5} = problem.equations{5}/model.scaling.qOs;
            problem.equations{6} = problem.equations{6}/model.scaling.qGs;
            problem.equations{7} = problem.equations{7}/model.scaling.bhp;
            
            
            
            eqJac = vertcat(problem.equations{:});
            Val = eqJac.val;
            Jac = cell2mat(eqJac.jac(1:3))*JacX0;
            
            varargout{1} = Val;
            varargout{2} = Jac;
end

    end
    function [varargout] = testGetEquationState(model, state0, stateVector,wellSol, dt, drivingForces,residualScale,varargin)
        
        
        if nargout == 1
            varargin = [varargin,{'resOnly',true}];
            
            [stateMrst] = model.toMRSTStates(stateVector);
            stateMrst.wellSol = wellSol;
            [problem] = getEquations(model, state0, stateMrst, dt, drivingForces,varargin{:},'reverseMode',false);
            
            problem.equations{1} = problem.equations{1}./residualScale.o;
            problem.equations{2} = problem.equations{2}./residualScale.w;
            problem.equations{3} = problem.equations{3}./residualScale.g;
            problem.equations{4} = problem.equations{4}/model.scaling.qWs;
            problem.equations{5} = problem.equations{5}/model.scaling.qOs;
            problem.equations{6} = problem.equations{6}/model.scaling.qGs;
            problem.equations{7} = problem.equations{7}/model.scaling.bhp;
            
            
            
            Val = vertcat(problem.equations{:});
            
            varargout{1} = Val;
            
        else
            varargin = [varargin,{'resOnly',false}];
            
            [stateMrst,JacX0] = model.toMRSTStates(stateVector);
            stateMrst.wellSol = wellSol;
            [problem] = getEquations(model, state0, stateMrst, dt, drivingForces,varargin{:},'reverseMode',false);
            
            
            problem.equations{1} = problem.equations{1}./residualScale.o;
            problem.equations{2} = problem.equations{2}./residualScale.w;
            problem.equations{3} = problem.equations{3}./residualScale.g;
            problem.equations{4} = problem.equations{4}/model.scaling.qWs;
            problem.equations{5} = problem.equations{5}/model.scaling.qOs;
            problem.equations{6} = problem.equations{6}/model.scaling.qGs;
            problem.equations{7} = problem.equations{7}/model.scaling.bhp;
            
            
            
            eqJac = vertcat(problem.equations{:});
            Val = eqJac.val;
            Jac = cell2mat(eqJac.jac(1:3))*JacX0;
            
            varargout{1} = Val;
            varargout{2} = Jac;
            
        end
        
    end
    
    function [varargout] = testGetEquationWell(model, state0, state,wellSolVector, dt, drivingForces,residualScale,varargin)
        
        
        if nargout == 1
            varargin = [varargin,{'resOnly',true}];
            
            [wellSol] = model.toWellSol(wellSolVector,drivingForces.Wells);
            state.wellSol = wellSol;
            [problem] = getEquations(model, state0, state, dt, drivingForces,varargin{:},'reverseMode',false);
            
            problem.equations{1} = problem.equations{1}./residualScale.o;
            problem.equations{2} = problem.equations{2}./residualScale.w;
            problem.equations{3} = problem.equations{3}./residualScale.g;
            problem.equations{4} = problem.equations{4}/model.scaling.qWs;
            problem.equations{5} = problem.equations{5}/model.scaling.qOs;
            problem.equations{6} = problem.equations{6}/model.scaling.qGs;
            problem.equations{7} = problem.equations{7}/model.scaling.bhp;
            
            
            Val = vertcat(problem.equations{:});
            
            varargout{1} = Val;
            
        else
            varargin = [varargin,{'resOnly',false}];
            
            [wellSol,JacW] = model.toWellSol(wellSolVector,drivingForces.Wells);
            state.wellSol = wellSol;
            [problem] = getEquations(model, state0, state, dt, drivingForces,varargin{:},'reverseMode',false);
            
            
            problem.equations{1} = problem.equations{1}./residualScale.o;
            problem.equations{2} = problem.equations{2}./residualScale.w;
            problem.equations{3} = problem.equations{3}./residualScale.g;
            problem.equations{4} = problem.equations{4}/model.scaling.qWs;
            problem.equations{5} = problem.equations{5}/model.scaling.qOs;
            problem.equations{6} = problem.equations{6}/model.scaling.qGs;
            problem.equations{7} = problem.equations{7}/model.scaling.bhp;
            
            
            eqJac = vertcat(problem.equations{:});
            Val = eqJac.val;
            Jac = cell2mat(eqJac.jac(4:7))*JacW;
            
            varargout{1} = Val;
            varargout{2} = Jac;
            
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
