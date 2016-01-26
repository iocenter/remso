function [ state ] = mrstStateConvexComb( alpha,state1,state2 )
% TODO: wells too!


withWells = (isfield(state1,'wellSol') || isfield(state2,'wellSol'));
if withWells && isfield(state1,'wellSol') && ~isfield(state2,'wellSol')
    state2.wellSol = state1.wellSol; %% assume both have same values
end


state = state2;


state.pressure = alpha * state1.pressure + (1-alpha)*state2.pressure;
state.s = alpha * state1.s + (1-alpha)*state2.s;
if isfield(state,'rs')
    state.rs = alpha * state1.rs + (1-alpha)*state2.rs;
end
if isfield(state,'rv')
    state.rv = alpha * state1.rv + (1-alpha)*state2.rv;
end
if isfield(state,'status')
    state = rmfield(state,'status');
end

if withWells
    for w = numel(state.wellSol)
        state.wellSol(w).qWs = alpha * state1.wellSol(w).qWs + (1-alpha)* state2.wellSol(w).qWs;
        state.wellSol(w).qOs = alpha * state1.wellSol(w).qOs + (1-alpha)* state2.wellSol(w).qOs;
        state.wellSol(w).qGs = alpha * state1.wellSol(w).qGs + (1-alpha)* state2.wellSol(w).qGs;
        state.wellSol(w).bhp = alpha * state1.wellSol(w).bhp + (1-alpha)* state2.wellSol(w).bhp;
    end
end


end

