function [ duC,dx,dv,xi,lowActive,upActive,mu,violationH,qpVAl,dxN,dvN,s] = qpStep(M,Bc,ldu,udu,ax,Ax,ldx,udx,av,Av,ldv,udv,varargin )
% Solves the Convex - QP problem:
%
%      qpVAl = min 1/2 duC'*M*du + Bc * duC
%           s,duC,xi st.
%                    ldx - s <= (1-xi)*ax + Ax * duC  <= udx + s
%                    ldv - s <= (1-xi)*av + Av * duC  <= udv + s
%                    ldu <= duC <= udu
%                    xi = xiM
%                    s = sM
%
% where xiM = min xi + bigM s
%          s,duC,xi st.
%                   ldx - s <= (1-xi)*ax + Ax * duC <= udx + s
%                   ldv - s <= (1-xi)*av + Av * duC <= udv + s
%                   lbu <= duC <= udu
%                   0 <= xi <= 1
%                   0 <= s  <= 1/bigM
%
%  The problem is solved iteratively.  We start with a subset of the
%  constraints. Given the solution of the iteration i, the violated
%  constratins are identified and included to the prevoius problem and re-solved.
%
% SYNOPSIS:
%  [ duC,dx,dv,lowActive,upActive,mu,s,violationH,qpVAl] = prsqpStep(M,Bc,u,lbu,ubu,Ax,ldx,udx,Av,ldv,udv )
%  [ duC,dx,dv,lowActive,upActive,mu,s,violationH,qpVAl] = prsqpStep(M,Bc,u,lbu,ubu,Ax,ldx,udx,Av,ldv,udv , 'pn', pv, ...)
% PARAMETERS:
%
%  Vectors and matrices to form the problem above, some of the represented
%  with cellarrays, related to the strucutre of the problem
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters. The
%             supported options are:
%
%   qpDebug - if true print information
%
%   lowActive, upActive - represent a guess for the set of tight constraints.
%
%   feasTol - Tolerance to judge the constraint activity
%             Judge the optimality condition quality and debug
%             Judge the feasibility of the problem.
%
%   ci - control incidence function
%
%   maxQpIt - Maximum number of QP iterations
%
%   it - remso it number for debug
%
% RETURNS:
%          listed the non-obvious (see the problem definition above)
%
%   dx = Ax * duC
%   dv = Av * duC
%   lowActive,upActive - tight constraints at the optimal solution.
%   mu - dual variables for the inequality constraints of the QP problem
%   violationH = Vector of norm1 of the constraint violation per iteration.
%
% SEE ALSO:
%
%


opt = struct('qpDebug',true,'lowActive',[],'upActive',[],'feasTol',1e-6,'ci',[],'maxQpIt',20,'it',0,'bigM',1e9);
opt = merge_options(opt, varargin{:});

if isempty(opt.ci)
    opt.ci = @(kk)controlIncidence([],kk);
end

withAlgs = ~isempty(ldv);


if opt.qpDebug
    if opt.it <= 1
        fid = fopen('logQP.txt','w');
        fidCplex = fopen('logCplex.txt','w');
    else
        fid = fopen('logQP.txt','a');
        fidCplex = fopen('logCplex.txt','a');
        
    end
    fprintf(fid,'********************* Iteration %3.d ********************\n',opt.it);
    fprintf(fid,'it l1-MinV AddCons LP-TIME ST l1-CurV NewCons QP-TIME ST\n');
    
    fprintf(fidCplex,'********************* Iteration %3.d ********************\n',opt.it);
    
    DisplayFunc = @(x) fprintf(fidCplex,[x,'\n']);
else
    DisplayFunc = [];
end


nu = numel(ldu{1});
nx = numel(udx{1});

if withAlgs
    nv = numel(udv{1});
end

nuH = sum(cellfun(@numel,ldu));

B =cell2mat(Bc);

du = [];
dx = [];
dv = [];

minViolation = 0;
violationH = [];


lowActiveP = opt.lowActive;
upActiveP = opt.upActive;

linesAdded.ub.x = upActiveP.x;
linesAdded.lb.x = lowActiveP.x;
if withAlgs
    linesAdded.ub.v = upActiveP.v;
    linesAdded.lb.v = lowActiveP.v;
end

newCons = cell(opt.maxQpIt+1,1);
nc = cell(opt.maxQpIt+1,1);

% mantain a record of the constraints added at each iteration
newCons{1} = linesAdded;

% count the number of new constraints
nc{1}.ub.x = sum(cellfun(@sum,linesAdded.ub.x));
nc{1}.lb.x = sum(cellfun(@sum,linesAdded.lb.x));
if withAlgs
    nc{1}.ub.v = sum(cellfun(@sum,linesAdded.ub.v));
    nc{1}.lb.v = sum(cellfun(@sum,linesAdded.lb.v));
end

violation = [];
solved = false;


nAddRows = 0;

for k = 1:opt.maxQpIt
    
    % extract the current constraints lines to add in this iteration
    [ Alow.x,blow.x,xRlow ] = buildActiveConstraints(Ax,ldx,newCons{k}.lb.x,-1,'R',ax);
    [ Aup.x,bup.x,xRup ]   = buildActiveConstraints(Ax,udx,newCons{k}.ub.x,1,'R',ax);
    if withAlgs
        [ Alow.v,blow.v,vRlow ] = buildActiveConstraints(Av,ldv,newCons{k}.lb.v,-1,'R',av);
        [ Aup.v,bup.v,vRup ]   = buildActiveConstraints(Av,udv,newCons{k}.ub.v,1,'R',av);
    end
    
    
    % Merge constraints in one matrix !
    if withAlgs
        Aq = cell2mat([[xRlow;xRup;vRlow;vRup],...
             [cellfun(@(xi)-ones(sum(xi),1),newCons{k}.lb.x,'UniformOutput',false);...
              cellfun(@(xi)-ones(sum(xi),1),newCons{k}.ub.x,'UniformOutput',false);...
              cellfun(@(xi)-ones(sum(xi),1),newCons{k}.lb.v,'UniformOutput',false);...
              cellfun(@(xi)-ones(sum(xi),1),newCons{k}.ub.v,'UniformOutput',false)],...
            [Alow.x;Aup.x;Alow.v;Aup.v]]);
        bq = cell2mat([blow.x;bup.x;blow.v;bup.v]);
    else
        Aq = cell2mat([[xRlow;xRup],...
             [cellfun(@(xi)-ones(sum(xi),1),newCons{k}.lb.x,'UniformOutput',false);...
              cellfun(@(xi)-ones(sum(xi),1),newCons{k}.ub.x,'UniformOutput',false)],...
            [Alow.x;Aup.x]]);
        bq = cell2mat([blow.x;bup.x]);
    end
    
    
    
    if k == 1
        P = Cplex('LP - QP');
        P.DisplayFunc = DisplayFunc;
        P.Param.qpmethod.Cur = 6;
        P.Param.lpmethod.Cur = 6;
        P.Param.emphasis.numerical.Cur = 1;
        
        
        % Add control variables
        P.addCols([zeros(nuH+2,1)], [], [0;0;cell2mat(ldu)], [1;1/opt.bigM;cell2mat(udu)]);  % objective will be set up again later
        P.Param.timelimit.Cur = 7200;
        P.addRows(0,[1,0,zeros(1,nuH)],1);
        P.addRows(0,[0,1,zeros(1,nuH)],1/opt.bigM);
    else
        P.Model.rhs(1) = 1;
        P.Model.lhs(1) = 0;
        P.Model.rhs(2) = 1/opt.bigM;
        P.Model.lhs(2) = 0;   
    end
    
    if opt.qpDebug
        fprintf(fidCplex,'************* LP %d  ******************\n',k);
    end
    
    
    % now add the new rows
    P.addRows(-inf(size(bq)),Aq,bq);
    
    
    % keep track of the number of constraints added so far
    nAddRows = nAddRows + numel(bq);
    
    
    % Set up the LP objective
    P.Model.Q = [];
    LPobj = [-1;opt.bigM;zeros(nuH,1)];
    P.Model.obj = LPobj;
    
    
    tic;
    P.solve();
    lpTime = toc;
    
    lpTime2 = 0;
    if (P.Solution.status ~= 1)
        method = P.Solution.method;
        status = P.Solution.status;
        if P.Solution.method == 2
            P.Param.lpmethod.Cur = 4;
        else
            P.Param.lpmethod.Cur = 2;
        end
        tic;
        P.solve();
        lpTime2 = toc;
        
        if (P.Solution.status ~= 1)
            if opt.qpDebug
                fprintf(fid,['Problems solving LP\n' ...
                    'Method ' num2str(method) ...
                    ' Status ' num2str(status) ...
                    ' Time %1.1e\n' ...
                    'Method ' num2str(P.Solution.method) ...
                    ' Status ' num2str(P.Solution.status)...
                    ' Time %1.1e\n'],lpTime,lpTime2) ;
            end
            warning(['Problems solving LP\n Method ' num2str(method) ' Status ' num2str(status) '\n Method ' num2str(P.Solution.method)  ' Status ' num2str(P.Solution.status)]);
        end
        
        P.Param.lpmethod.Cur = 6;
    end
    
    % Determine the value of 'xibar', for the current iteration ,see the problem difinition above
    xibar = P.Solution.x(1);
    sM = P.Solution.x(2);
    
    if opt.qpDebug
        fprintf(fid,'%2.d %1.1e %1.1e %1.1e %2.d ',k,1-xibar,nAddRows,lpTime+lpTime2,P.Solution.status) ;
        fprintf(fidCplex,'************* QP %d  *******************\n',k);
    end
    
    % set up the qp objective
    P.Model.Q = blkdiag(1,1,M);
    P.Model.obj = [-1;opt.bigM;B'];
    P.Model.rhs(1) = xibar;
    P.Model.lhs(1) = xibar;
    P.Model.rhs(2) = sM;
    P.Model.lhs(2) = sM;
    
    tic;
    P.solve();
    QpTime = toc;
    
    QpTime2 = 0;
    if (P.Solution.status ~= 1)
        method = P.Solution.method;
        status = P.Solution.status;
        if P.Solution.method == 2
            P.Param.qpmethod.Cur = 4;
        else
            P.Param.qpmethod.Cur = 2;
        end
        tic;
        P.solve();
        QpTime2 = toc;
        if (P.Solution.status ~= 1)
            if opt.qpDebug
                fprintf(fid,['Problems solving QP\n' ...
                    'Method ' num2str(method) ...
                    ' Status ' num2str(status) ...
                    ' Time %1.1e\n' ...
                    'Method ' num2str(P.Solution.method) ...
                    ' Status ' num2str(P.Solution.status)...
                    ' Time %1.1e\n'],QpTime,QpTime2) ;
            end
            warning(['\nProblems solving QP\n Method ' num2str(method) ' Status ' num2str(status) '\n Method ' num2str(P.Solution.method)  ' Status ' num2str(P.Solution.status)]);
        end
        
        P.Param.qpmethod.Cur = 6;
    end
    
    
    
    du = P.Solution.x(3:end);
    
    
    duC = toStructuredCells(du,nu);
    
    dxN = cellmtimes(Ax,duC,'lowerTriangular',true,'ci',opt.ci);
    if withAlgs
        dvN = cellmtimes(Av,duC,'lowerTriangular',true,'ci',opt.ci);
    else
        dvN = [];
    end
    
    du = duC;
    dx = cellfun(@(z,dz)xibar*z+dz,ax,dxN,'UniformOutput',false);
    if withAlgs
        dv = cellfun(@(z,dz)xibar*z+dz,av,dvN,'UniformOutput',false);
    end
    
    % Check which other constraints are infeasible
    [feasible.x,lowActive.x,upActive.x,violation.x ] = checkConstraintFeasibility(dx,ldx,udx,'primalFeasTol',opt.feasTol ) ;
    if withAlgs
        [feasible.v,lowActive.v,upActive.v,violation.v ] = checkConstraintFeasibility(dv,ldv,udv,'primalFeasTol',opt.feasTol  );
    end
    
    
    % debugging purpouse:  see if the violation is decreasing!
    ineqViolation = violation.x;
    if withAlgs
        ineqViolation = ineqViolation + violation.v;
    end
    
    violationH = [violationH,ineqViolation];
    
    
    
    % Determine the new set of contraints to be added
    [newCons{k+1}.lb.x,nc{k+1}.lb.x] =  newConstraints(linesAdded.lb.x, lowActive.x);
    [newCons{k+1}.ub.x,nc{k+1}.ub.x] =  newConstraints(linesAdded.ub.x, upActive.x);
    if withAlgs
        [newCons{k+1}.lb.v,nc{k+1}.lb.v] =  newConstraints(linesAdded.lb.v, lowActive.v);
        [newCons{k+1}.ub.v,nc{k+1}.ub.v] =  newConstraints(linesAdded.ub.v, upActive.v);
    end
    
    % determine how many new constraints are being added in total
    newC = 0;
    bN = fieldnames(nc{k+1});
    for bi = 1:numel(bN)
        vN = fieldnames(nc{k+1}.(bN{bi}));
        for vi = 1:numel(vN);
            newC  = newC + nc{k+1}.(bN{bi}).(vN{vi});
        end
    end
    
    
    if opt.qpDebug
        fprintf(fid,'%1.1e %1.1e %1.1e %2.d\n',ineqViolation,newC,QpTime+QpTime2,P.Solution.status) ;
    end
    
    % if we cannot add more constraints, so the problem is solved!
    if newC == 0
        if minViolation > opt.feasTol
            if opt.qpDebug
                fprintf(fid,'Irreductible constraint violation l1 norm: %e \n',minViolation) ;
            end
        end
        solved = true;
        break;
    end
    
    
    linesAdded.lb.x = unionActiveSets(linesAdded.lb.x,lowActive.x);
    linesAdded.ub.x = unionActiveSets(linesAdded.ub.x,upActive.x);
    if withAlgs
        linesAdded.lb.v = unionActiveSets(linesAdded.lb.v,lowActive.v);
        linesAdded.ub.v = unionActiveSets(linesAdded.ub.v,upActive.v);
    end
    
end
if ~solved
    warning('Qps were not solved within the iteration limit')
end

xi = 1-xibar;
s = sM;

% extract the dual variables:
to = 2;  %% skip xibar and sM dual constraint!
for j = 1:k
    from = to + 1;
    to = from + nc{j}.lb.x-1;
    [r] = extractCompressIneq(-P.Solution.dual(from:to),newCons{j}.lb.x,nx);
    if j==1
        mu.lb.x = r;
    else
        mu.lb.x = cellfun(@(x1,x2)x1+x2,mu.lb.x,r,'UniformOutput',false);
    end
    
    from = to + 1;
    to = from + nc{j}.ub.x-1;
    [r] = extractCompressIneq(-P.Solution.dual(from:to),newCons{j}.ub.x,nx);
    if j==1
        mu.ub.x = r;
    else
        mu.ub.x = cellfun(@(x1,x2)x1+x2,mu.ub.x,r,'UniformOutput',false);
    end
    
    if withAlgs
        from = to + 1;
        to = from + nc{j}.lb.v-1;
        [r] = extractCompressIneq(-P.Solution.dual(from:to),newCons{j}.lb.v,nv);
        if j==1
            mu.lb.v = r;
        else
            mu.lb.v = cellfun(@(x1,x2)x1+x2,mu.lb.v,r,'UniformOutput',false);
        end
        
        from = to + 1;
        to = from + nc{j}.ub.v-1;
        [r] = extractCompressIneq(-P.Solution.dual(from:to),newCons{j}.ub.v,nv);
        if j==1
            mu.ub.v = r;
        else
            mu.ub.v = cellfun(@(x1,x2)x1+x2,mu.ub.v,r,'UniformOutput',false);
        end
    end
end
if to ~= size(P.Model.rhs,1)
    error('checkSizes!');
end

%%%  First order optimality
%norm(P.Model.Q * P.Solution.x + P.Model.obj -(P.Model.A)' * P.Solution.dual - P.Solution.reducedcost)


% extract dual variables with respect to the controls
mu.ub.u = toStructuredCells(max(-P.Solution.reducedcost(3:nuH+2),0),nu);
mu.lb.u = toStructuredCells(max( P.Solution.reducedcost(3:nuH+2),0),nu);

if opt.qpDebug
    if withAlgs
        optCheck = cell2mat(duC)'*M + cell2mat(Bc) + ...
            cell2mat(cellmtimesT( cellfun(@(x1,x2)x1-x2,mu.ub.x,mu.lb.x,'UniformOutput',false),Ax ,'lowerTriangular',true,'ci',opt.ci)) + ...
            cell2mat(cellmtimesT( cellfun(@(x1,x2)x1-x2,mu.ub.v,mu.lb.v,'UniformOutput',false),Av,'lowerTriangular',true,'ci',opt.ci)) + ...
            cell2mat(             cellfun(@(x1,x2)x1-x2,mu.ub.u,mu.lb.u,'UniformOutput',false))';
    else
        optCheck = cell2mat(duC)'*M + cell2mat(Bc) + ...
            cell2mat(cellmtimesT( cellfun(@(x1,x2)x1-x2,mu.ub.x,mu.lb.x,'UniformOutput',false),Ax ,'lowerTriangular',true,'ci',opt.ci)) + ...
            cell2mat(             cellfun(@(x1,x2)x1-x2,mu.ub.u,mu.lb.u,'UniformOutput',false))';
    end
    
    optNorm = norm(optCheck);
    fprintf(fid,'Optimality norm: %e \n',optNorm) ;
    if optNorm > opt.feasTol*10
        warning('QP optimality norm might be to high');
    end
end


qpVAl = P.Solution.objval;

if opt.qpDebug
    
    
    fclose(fid);
end

end



