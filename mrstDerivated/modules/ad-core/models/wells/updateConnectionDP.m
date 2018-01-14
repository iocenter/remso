function sol = updateConnectionDP(wellmodel, model, sol)
% Explicit update of hydrostatic pressure difference between bottom hole
% and connections based on phase distrubution alonw well-bore.
%
% SYNOPSIS:
%   sol = updateConnectionDP(wellmodel, model, sol)
%
% PARAMETERS:
%   wellmodel   - Simulation well model.
%   model       - Simulation model.
%   sol         - List of current well solution structures
%
% RETURNS:
%   sol         - Well solution structures with updated field 'cdp'
%
% SEE ALSO:
%   WellModel, computeWellContributionsNew.

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
% Explicit update of hydrostatic pressure difference between bottom hole
% and connections.
% input:
% sol : well-solutions with
%{
Changes from https://github.com/iocenter/remso/ 4f8fa54a92c14b117445ad54e9e5dd3a0e47a7f5
  Make operations compatible with ADI objects

%}
W = wellmodel.W;
b = wellmodel.bfactors;
rhos = wellmodel.surfaceDensities;
rMax = wellmodel.maxComponents;

if ~iscell(sol(1).cqs)
    sol = arrayfun(@(si) subsasgn(si,struct('type',{'.'},'subs',{'cqs'}),... %sol(it).cqs =
                         mat2cell(si.cqs,size(si.cqs,1),ones(1,size(si.cqs,2)))),...  %mat2cell(sol(it).cqs)
                         sol);
end

nConn       = cellfun(@numel, {W.cells})'; % # connections of each well
perf2well   = rldecode((1:numel(W))', nConn);

numPh = numel(b);

% if nargin < 6
%     model = getModel(size(b,2), size(r,2));
% end
% actPh = getActivePhases(model);
%[isActive, actPh] = model.getActivePhases();


for k = 1:numel(sol);
    s  = sol(k);
    w  = W(k);

    if ~isfield(w, 'topo')
        nperf = numel(w.cells);
        w.topo = [(0:(nperf-1))', (1:nperf)'];
    end

    qs = s.cqs; %volumetric in-flux at standard conds
    perfInx = (perf2well == k);
    bk      = cellfun(@(bi)bi(perfInx),b,'UniformOutput',false);
    if ~isempty(rMax)
        rkMax   = cellfun(@(rMaxi)rMaxi(perfInx),rMax,'UniformOutput',false);
    else
        rkMax   = [];
    end

    C = wb2in(w);            % mapping wb-flux to in-flux
    wbqs  = cellfun(@(qsi)abs(C\qsi),qs,'UniformOutput',false);       % solve to get well-bore fluxes at surface conds
    wbqst = repmat(speye(numel(double(wbqs{1}))),1,numel(wbqs)  )*vertcat(wbqs{:});  % wbqst = sum(wbqs, 2);   % total wb-flux at std conds
    % if flux is zero - just use compi
    zi = double(wbqst) == 0;
    if any( zi )
        wbqsZ = cellfun(@(ci)ones(nnz(zi),1)*ci,...                       %wbqsZ =  ones(nnz(zi),1)*w.compi
                                  num2cell(w.compi(1:numel(wbqs))),...
                                  'UniformOutput',false);
        wbqs = cellfun(@(wbqsi,wbqsZi)subsasgn(wbqsi,struct('type','()','subs',{{zi}}),wbqsZi),wbqs,wbqsZ,'UniformOutput',false);      % wbqs(zi) = wbqsZ                  
        wbqst(zi) = repmat(speye(nnz(zi)),1,numel(wbqsZ)  )*vertcat(wbqsZ{:});  %wbqst(zi,:) = sum(wbqsZ, 2);
    end
    % Compute mixture at std conds:
    mixs = cellfun(@(wbqsi)wbqsi./wbqst,wbqs,'UniformOutput',false);
    % compute volume ratio Vr/Vs
    volRat = compVolRat(mixs, bk, rkMax, model);
    % Mixture density at connection conds (by using static b's)
    rhoMix = cell2mat(arrayfun(@(ri)ri*speye(numel(double(mixs{1}))),rhos,'UniformOutput',false))*vertcat(mixs{:}) ./volRat;%(mixs*rhos(:))./volRat;
    % rhoMix is now density between neighboring segments given by
    % topo(:,1)-topo(:,2) computed by using conditions in well-cell
    % topo(:,2). This is probably sufficiently accurate.

    % get dz between segment nodes and bh-node1
    dpt = [0; w.dZ];
    dz  = diff(dpt);
    g   = norm(gravity);
    ddp = g*rhoMix.*dz; % p-diff between connection neighbors
    % well topology assumes we can traverse from top down, but add a loop
    % just in case crazy ordering.
    cdp    = ddp*nan;
    cdp(1) = ddp(1);
    its = 0; maxIts = 100;
    while and(any(isnan(cdp)), its<maxIts)
        its = its +1;
        for cnr = 2:numel(double(cdp))
            cdp(w.topo(cnr,2)) = cdp(w.topo(cnr,1)) + ddp(cnr);
        end
    end
    if its == maxIts
        error(['Problem with topology for well: ', s.name, '. Segments appear not to be connected'])
    end
    sol(k).cdp = cdp;
end
end


function C = wb2in(w)
    conn = w.topo(2:end, :);
    % Number of connections between perforations
    nconn = size(conn, 1);
    % Number of perforations
    nperf = numel(w.cells);
    
    if nconn + 1 ~= nperf
        warning(['Mismatch between connection count (', num2str(nconn+1),...
                ') and perforation count (', num2str(nperf), '). Well model', ...
                'Does not appear to be a tree.']);
    end

    id = (1:nperf)';
    % First identity, then honor topology.
    ii = [id; conn(:, 1)];
    jj = [id; conn(:, 2)];
    vv = [ones(nperf, 1); -ones(nconn, 1)]; 

    C = sparse(ii, jj, vv, nperf, nperf);
end

function volRat = compVolRat(mixs, b, rMax, model)
%
x = mixs;
dg = isprop(model, 'disgas') && model.disgas;
vo = isprop(model, 'vapoil') && model.vapoil;

if dg || vo
    [~, isgas] = model.getVariableField('sg');
    [~, isoil] = model.getVariableField('so');
    
    both = find(isgas | isoil);
    
    g = mixs{isgas};
    o = mixs{isoil};
    if isa(model, 'ThreePhaseBlackOilModel')
        % Vapoil/disgas
        gor = abs(g./o);
        gor(or(isnan(gor),isinf(gor))) = inf;
        rs = min(rMax{1}, gor);
        ogr = abs(o./g);
        ogr(or(isnan(ogr),isinf(ogr))) = inf;
        rv = min(rMax{2}, ogr);
        d = 1-rs.*rv;
        x{isgas} = (x{isgas} - rs.*o)./d;
        x{isoil} = (x{isoil} - rv.*g)./d;
        x{both} = x{both}.*(x{both}>0);
    else
        % Only gas dissolution
        x{isgas} = x{isgas} - rMax{1}.*o;
        x{isgas} = x{isgas}.*(x{isgas}>0);
    end
end

ratio = cellfun(@(xi,bi)xi./bi,x,b,'UniformOutput',false);
volRat = repmat(speye(numel(double(ratio{1}))),1,numel(ratio)  )*vertcat(ratio{:});  %  sum(ratio ,2);
end



