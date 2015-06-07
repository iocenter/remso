function [grad] = runForwardGradient(G, S, W, rock, fluid, simRes,schedule, controls,sRHS,uRHS,obj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%{

P(s0,x,u) = A(s0) x - b(u) = 0   %% PressureSystem
T(s,s0,x) = s-s0-invPV*(A(x)fw(s)+q(v)) = 0

% b may depend on the gravity, but it is disregarded by the adjoint method.
I'm not sure if it is safe to disregard it, but to be sure, use this method
with 'gravity off' and no capillary pressure to be sure that the gradients
are correct!

%}

bhpWells  = strcmpi('bhp' , { W.type });
rateWells = strcmpi('rate', { W.type });
WS         = [ W.S ];
Dw        = blkdiag( WS.D );
DwD       = Dw(:, bhpWells);
[A_N, ~, A_D, ~] = controls2Wells(W, schedule, controls);
numW = numel(W);

Bqtot = sparse(numW,numW);

Bqtot(rateWells,:) = A_N{1};
Bpw = DwD*A_D{1};


numCF    = size(G.cells.faces, 1);
numC     = G.cells.num;
numF     = G.faces.num;
numU     = numel(controls.well);

sizeSeeds = size(uRHS,2);
BU = -[sparse(numCF,sizeSeeds);Bpw*uRHS;sparse(numC,sizeSeeds);sparse(numF,sizeSeeds);Bqtot*uRHS];

% Assuming a single step
curStep = 2;

dt      = simRes(curStep).timeInterval * [-1 1]';

% Generate RHS, that is f-part (rest is zero)
numC    = G.cells.num;
PV      = G.cells.volumes.*rock.poro;
invPV  = 1./PV;
invDPV  = spdiags(invPV, 0, numC, numC);
%f_w     = fluid.krw(simRes(curStep).resSol) ./
%fluid.Lt(simRes(curStep).resSol); % fractinal flow, f(s^n)
%f_w     = fluid.fw( simRes(curStep).resSol);

%{
mu  = fluid.properties(simRes(curStep).resSol);
sat = fluid.saturation(simRes(curStep).resSol);
kr  = fluid.relperm(sat, simRes(curStep).resSol);
mob = bsxfun(@rdivide, kr, mu);
f_w = mob(:,1) ./ sum(mob,2);
%}
[mob, dmob] = mobilities(simRes(curStep-1).resSol, fluid);
Lt          = sum(mob, 2);

if strcmp(S.type, 'hybrid')
    solver = 'hybrid';
else
    solver = 'mixed';
end

% Generate right-hand-side
%dim     = fluid.dLtInv(simRes(curStep).resSol);
dim     = -sum(dmob, 2) ./ Lt.^2;   % d/ds (1/Lt)


cellNo  = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
S.C     = sparse(1:numel(cellNo), cellNo, 1);
sgn     = 2*double(G.faces.neighbors(G.cells.faces(:,1),1) == cellNo) - 1;
v       = sgn .* simRes(curStep).resSol.flux(G.cells.faces(:,1));   % v^{n+1}
%l_v     = bsxfun(@times,sgn, adjRes(curStep+1).resSol.flux(G.cells.faces(:,1),:));   % lam_v^{n+1}
vFs0     = -((S.C*sRHS)'*spdiags( (S.C*dim).*v , 0, numCF, numCF)*S.B)';

nperf = cellfun(@numel, { W.cells });
snperf = sum(cellfun(@numel, { W.cells }));
perfI = rldecode(1:numW, nperf, 2);
qFs0 = sparse(snperf,sizeSeeds);
for wellNr = 1:numW
    w   = W(wellNr);
    q   = simRes(curStep).resSol.wellSol(wellNr).flux; % q_w^{n+1}
    qFs0(perfI==wellNr,:) = ((w.S.C*sRHS)'*bsxfun(@times, ( (w.S.C*dim).*q ),(w.S.B) ))';
end


Bs0 = [vFs0;qFs0;sparse(numC,sizeSeeds);sparse(numF+numU,sizeSeeds)];


b = BU+Bs0;
b = mat2cell(b,[numCF+snperf,numC,numF+numU],sizeSeeds)';


% Solve linear system based on s^{n-1}
resSolFwd = solveIncompFlowLocal(simRes(curStep-1).resSol, G, S, fluid, ...
    'gravityOff',true,...
    'wells', W, 'rhs', b, 'Solver', solver);

%{

% Solve linear system based on s^{n-1}
resSolFwd2 = solveIncompFlowLocal(simRes(curStep-1).resSol, G, S, fluid, ...
    'gravityOff',true,...
    'wells', W, 'rhsS', sRHS, 'Solver', solver);


norm(resSolFwd2.pressure - resSolFwd.pressure,inf)
norm(resSolFwd2.flux - resSolFwd.flux,inf)
norm(vertcat(resSolFwd2.wellSol.flux) - vertcat(resSolFwd.wellSol.flux),inf)
norm(vertcat(resSolFwd2.wellSol.pressure) - vertcat(resSolFwd.wellSol.pressure),inf)

%}


%% ++++++++++++++++++  Transport part
[mob, dmob] = mobilities(simRes(curStep).resSol, fluid);
Lt          = sum(mob, 2);
f           = bsxfun(@rdivide, mob, Lt);
f_w         = f(:,1);

Dfw_all     = (f(:,2).*dmob(:,1) - f(:,1).*dmob(:,2)) ./ Lt;   % Chain rule.
DDf         = spdiags(Dfw_all, 0, numC, numC);

At      = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol,simRes(curStep).wellSol, 'Transpose', true);

% this matrix shoudln't be symmetric in general  --> due to (at least) At
systMat = speye(numC, numC) - dt * ( DDf * At * invDPV)';   % system matrix

dvF    = bsxfun(@times,sgn, resSolFwd.flux(G.cells.faces(:,1),:));   % lam_v^{n+1}

% Flux-matrix: A.i, A.j, A.qMinus
[A, qPluss, signQ] = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
    simRes(curStep).wellSol, 'VectorOutput', true);

dQPluss =  double( signQ > 0 );
dQMinus = -double( signQ < 0 );

dvFf    =  - dt*( bsxfun(@times,f_w(A.j),dvF)'*( S.C*(invDPV)) )'  ...
            + dt*( bsxfun(@times,(-f_w.*dQMinus + dQPluss).*invPV  ,S.C'*dvF ));

g_s = systMat\(     dvFf + sRHS      );
g_qw = vertcat(resSolFwd.wellSol(:).flux);


assert(obj.partials(2).v == 0, 'Assuming the target does not depend on v' );
grad = obj.partials(2).s* g_s+obj.partials(2).q_w *g_qw   ;

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
