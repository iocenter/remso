function [ duC,dx,dv,lowActive,upActive,mu,s,violation,qpVAl] = prsqpStep(M,Bc,u,lbu,ubu,Ax,ldx,udx,Av,ldv,udv,varargin )
% Solves the Convex - QP problem:
%
%      qpVAl = min 1/2 duC'*M*du + Bc * duC
%             duC,s st.
%                    sum(s.x)+sum(s.v) <= xi
%                    ldx - s.x <= Ax * duC <= udx + s.x
%                    ldv - s.v <= Av * duC <= udv + s.v
%                    lbu <= u + duC <= ubv
%                    s >= 0
%
% where xi = min sum(s.x)+sum(s.v)
%            duC,s st.
%                   ldx - s.x <= Ax * duC <= udx + s.x
%                   ldv - s.v <= Av * duC <= udv + s.v
%                   lbu <= u + duC <= ubv
%                   s >= 0
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


opt = struct('qpDebug',true,'lowActive',[],'upActive',[],'feasTol',1e-6,'ci',[],'maxQpIt',20,'it',0,'withAlgs',false);
opt = merge_options(opt, varargin{:});

if isempty(opt.ci)
    opt.ci = @(kk)controlIncidence([],kk);
end

withAlgs = opt.withAlgs;


if opt.qpDebug
    if opt.it <= 1
        fid = fopen('logQP.txt','w');
        fidCplex = fopen('logCplex.txt','w');
    else
        fid = fopen('logQP.txt','a');
        fidCplex = fopen('logCplex.txt','a');
        
    end
    fprintf(fid,'********************* Iteration %3.d ********************\n',opt.it);
    fprintf(fid,'it l1-MinV AddCons LP-TIME ST infCurV NewCons QP-TIME ST\n');
    
    fprintf(fidCplex,'********************* Iteration %3.d ********************\n',opt.it);
    
    DisplayFunc = @(x) fprintf(fidCplex,[x,'\n']);
else
    DisplayFunc = [];
end


uDims = cellfun(@numel,u);
xDims = cellfun(@numel,udx);

if withAlgs
    vDims = cellfun(@numel,udv);
end

uV = cell2mat(u);
nuH = numel(uV);

B =cell2mat(Bc);
s = [];

du = [];
dx = [];
dv = [];

minViolation = 0;
violationH = [];

elasticMode = true;

if isempty(ubu)
    udu = inf(nuH,1);
else
    udu = ubu-uV;
end
if isempty(lbu)
    ldu = -inf(nuH,1);
else
    ldu = lbu-uV;
end

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


nAddSlacks = 0;

for k = 1:opt.maxQpIt
    
    % extract the current constraints lines to add in this iteration
    [ Alow.x,blow.x ] = buildActiveConstraints(Ax,ldx,newCons{k}.lb.x,-1);
    [ Aup.x,bup.x ]   = buildActiveConstraints(Ax,udx,newCons{k}.ub.x,1);
    if withAlgs
        [ Alow.v,blow.v ] = buildActiveConstraints(Av,ldv,newCons{k}.lb.v,-1);
        [ Aup.v,bup.v ]   = buildActiveConstraints(Av,udv,newCons{k}.ub.v,1);
    end
    
    
    % Merge constraints in one matrix !
    if withAlgs
        Aq = cell2mat([Alow.x;Aup.x;Alow.v;Aup.v]);
        bq = cell2mat([blow.x;bup.x;blow.v;bup.v]);
    else
        Aq = cell2mat([Alow.x;Aup.x]);
        bq = cell2mat([blow.x;bup.x]);
    end
    
    if ~elasticMode
        error('not implemented')
        
    end
    if elasticMode
        
        if k == 1
            P = Cplex('LP - QP');
            P.DisplayFunc = DisplayFunc;
            P.Param.qpmethod.Cur = 6;
            P.Param.lpmethod.Cur = 6;
            P.Param.emphasis.numerical.Cur = 1;
            
            
            % Add control variables
            P.addCols(zeros(nuH,1), [], ldu, udu);  % objective will be set up again later
            P.Param.timelimit.Cur = 7200;
        else
            % remove the last constraint, sum(s) <= xi, this will be
            % reincorporated later with the new xi value
            P.delRows(numel(P.Model.rhs));
        end
        
        if opt.qpDebug
            fprintf(fidCplex,'************* LP %d  ******************\n',k);
        end
        %%% build LP
        nNewSlacks = size(bq,1);
        
        
        % add new variables  %%
        P.addCols(zeros(nNewSlacks,1), sparse(nAddSlacks,nNewSlacks), zeros(nNewSlacks,1), inf(nNewSlacks,1));
        
        
        % now add the new rows
        P.addRows(-inf(size(bq)),[Aq,sparse(nNewSlacks,nAddSlacks),-speye(nNewSlacks)],bq);
        
        
        % keep track of the number of constraints added so far
        nAddSlacks = nAddSlacks + nNewSlacks;
        
        
        % Set up the LP objective
        P.Model.Q = [];
        LPobj = [zeros(nuH,1);ones(nAddSlacks,1)];
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
        
        % Determine the value of 'xi', for the current iteration ,see the problem difinition above
        minViolation = max(sum([max(P.Model.A * P.Solution.x - P.Model.rhs,0);0]),...
            P.Solution.objval);
        
        if opt.qpDebug
            fprintf(fid,'%2.d %1.1e %1.1e %1.1e %2.d ',k,minViolation,nAddSlacks,lpTime+lpTime2,P.Solution.status) ;
            fprintf(fidCplex,'************* QP %d  *******************\n',k);
        end
        
        % set up the qp objective
        P.Model.Q = blkdiag(M,sparse(nAddSlacks,nAddSlacks));
        P.Model.obj = [B';zeros(nAddSlacks,1)];
        P.addRows(-inf,LPobj',minViolation);
        
        
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
        
        
        
        du = P.Solution.x(1:nuH);
        sSol = P.Solution.x(nuH+1:end);
        
    end
    
    duC = mat2cell(du,uDims,1);
    
    dx = cellmtimes(Ax,duC,'lowerTriangular',true,'ci',opt.ci);
    if withAlgs
        dv = cellmtimes(Av,duC,'lowerTriangular',true,'ci',opt.ci);
    else
        dv = [];
    end
    
    % Check which other constraints are infeasible
    [feasible.x,lowActive.x,upActive.x,violation.x ] = checkConstraintFeasibility(dx,ldx,udx,'primalFeasTol',opt.feasTol ) ;
    if withAlgs
        [feasible.v,lowActive.v,upActive.v,violation.v ] = checkConstraintFeasibility(dv,ldv,udv,'primalFeasTol',opt.feasTol  );
    end
    
    
    % debugging purpouse:  see if the violation is decreasing!
    ineqViolation = violation.x;
    if withAlgs
        ineqViolation = max(ineqViolation,violation.v);
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
                fprintf(fid,'Irreductible constraint violation norm_inf: %e \n',ineqViolation) ;
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

% extract the optimal slack variables
if elasticMode
    s.x = cellfun(@(x)zeros(size(x)),ldx,'UniformOutput',false);
    if withAlgs
        s.v = cellfun(@(x)zeros(size(x)),ldv,'UniformOutput',false);
    end
    
    to = 0;
    for j = 1:k
        
        from = to +1;
        to = from + nc{j}.lb.x-1;
        [sxl] = extractCompressIneq(sSol(from:to),newCons{j}.lb.x,xDims);
        
        from = to + 1;
        to = from + nc{j}.ub.x-1;
        [sxu] = extractCompressIneq(sSol(from:to),newCons{j}.ub.x,xDims);
        
        % sum here works just as simple assignment, since the variables
        % appears only once!
        if j == 1
            s.x = cellfun(@(x1,x2)x1+x2,sxl,sxu,'UniformOutput',false);
        else
            s.x = cellfun(@(x1,x2,x3)x1+x2+x3,s.x,sxl,sxu,'UniformOutput',false);
        end
        
        if withAlgs
            from = to +1;
            to = from + nc{j}.lb.v-1;
            [svl] = extractCompressIneq(sSol(from:to),newCons{j}.lb.v,vDims);
            
            from = to + 1;
            to = from + nc{j}.ub.v-1;
            [svu] = extractCompressIneq(sSol(from:to),newCons{j}.ub.v,vDims);
            
            if j == 1
                s.v = cellfun(@(x1,x2)x1+x2,svl,svu,'UniformOutput',false);
            else
                s.v = cellfun(@(x1,x2,x3)x1+x2+x3,s.v,svl,svu,'UniformOutput',false);
            end
        end
    end
    if to ~= size(sSol,1)
        error('checkSizes!');
    end
end

% extract the dual variables:
to = 0;
for j = 1:k
    from = to + 1;
    to = from + nc{j}.lb.x-1;
    [r] = extractCompressIneq(-P.Solution.dual(from:to),newCons{j}.lb.x,xDims);
    if j==1
        mu.lb.x = r;
    else
        mu.lb.x = cellfun(@(x1,x2)x1+x2,mu.lb.x,r,'UniformOutput',false);
    end
    
    from = to + 1;
    to = from + nc{j}.ub.x-1;
    [r] = extractCompressIneq(-P.Solution.dual(from:to),newCons{j}.ub.x,xDims);
    if j==1
        mu.ub.x = r;
    else
        mu.ub.x = cellfun(@(x1,x2)x1+x2,mu.ub.x,r,'UniformOutput',false);
    end
    
    if withAlgs
        from = to + 1;
        to = from + nc{j}.lb.v-1;
        [r] = extractCompressIneq(-P.Solution.dual(from:to),newCons{j}.lb.v,vDims);
        if j==1
            mu.lb.v = r;
        else
            mu.lb.v = cellfun(@(x1,x2)x1+x2,mu.lb.v,r,'UniformOutput',false);
        end
        
        from = to + 1;
        to = from + nc{j}.ub.v-1;
        [r] = extractCompressIneq(-P.Solution.dual(from:to),newCons{j}.ub.v,vDims);
        if j==1
            mu.ub.v = r;
        else
            mu.ub.v = cellfun(@(x1,x2)x1+x2,mu.ub.v,r,'UniformOutput',false);
        end
    end
end
to = to + 1;  %% skip least-feasibility constraint dual!
if to ~= size(P.Model.rhs,1)
    error('checkSizes!');
end

%%%  First order optimality
%norm(P.Model.Q * P.Solution.x + P.Model.obj -(P.Model.A)' * P.Solution.dual - P.Solution.reducedcost)


% extract dual variables with respect to the controls
mu.ub.u = mat2cell(max(-P.Solution.reducedcost(1:nuH),0),uDims,1);
mu.lb.u = mat2cell(max(P.Solution.reducedcost(1:nuH),0),uDims,1);

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



