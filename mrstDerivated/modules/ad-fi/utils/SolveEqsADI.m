function dx = SolveEqsADI(eqs, phi)

%
% Modification by Codas: 
% Deal with multiple RHS
%

useBasis = ~isempty(phi);

if useBasis
    % We have been provided a pod basis, use it to reduce equation sets
    % before solving.
    for i = 1:numel(eqs)
        eqs{i}.val = phi.basis{i}'*eqs{i}.val;
    end
end

numVars = cellfun(@(x) size(x.val,1), eqs)';
cumVars = cumsum(numVars); 
ii = [[1;cumVars(1:end-1)+1], cumVars];

eqs = cat(eqs{:});
tic

% Above CAT means '.jac' is a single element cell array.  Extract contents.
J = -eqs.jac{1};

if useBasis
    blkphi = blkdiag(phi.basis{:});
    J = blkphi'*J*blkphi;
end

[L,U,P,Q] = lu(J);
tmp = Q * (U \ (L \ (P * eqs.val))) ;

assert(all(all(isfinite(tmp))), 'Linear solver produced non-finite values!')

eqn = size(ii,1);
dx = cell(eqn,1);
for i = 1:eqn
    dx{i} = tmp(ii(i,1):ii(i,2),:);
end

if useBasis
    for i = 1:numel(dx)
        dx{i} = phi.basis{i}*dx{i};
    end
end
