function [ duC,dx,dv,ds,xi,lowActive,upActive,mu,violation,qpVAl,dxN,dvN,dsN,slack,k] = qpStep_R(M,g,w,ldu,udu,Aact1,predictor,constraintBuilder,ax,Ax,ldx,udx,av,Av,ldv,udv,as,As,lds,uds,varargin )
% Solves the Convex - QP problem:
%
%      qpVAl = min 1/2 duC'*M*du + (g + (1-xi*) w ) * duC
%             s,duC  st.
%                    ldx - s <= (1-xi*)*ax + Ax * duC  <= udx + s
%                    ldv - s <= (1-xi*)*av + Av * duC  <= udv + s
%                    ldu <= duC <= udu
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
%  [ duC,dx,dv,lowActive,upActive,mu,s,violationH,qpVAl] = prsqpStep(M,g,w,u,lbu,ubu,Ax,ldx,udx,Av,ldv,udv )
%  [ duC,dx,dv,lowActive,upActive,mu,s,violationH,qpVAl] = prsqpStep(M,g,w,u,lbu,ubu,Ax,ldx,udx,Av,ldv,udv , 'pn', pv, ...)
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


opt = struct('qpDebug',true,'lowActive',[],'upActive',[],'feasTol',1e-6,'ss',[],'maxQpIt',20,'it',0,'bigM',1e9,'condense',true,'lagFunc',[],'testQP',false);
opt = merge_options(opt, varargin{:});


%assert(opt.feasTol*10>1/opt.bigM,'Make sure bigM is big enough compared to feasTol');
opt.bigM = max(opt.bigM,10/opt.feasTol);

withAlgs = true;

xibar = 1;

if opt.qpDebug
    if opt.it <= 1
        fid = fopen('logQP.txt','w');
        fidCplex = fopen('logCplex.txt','w');
    else
        fid = fopen('logQP.txt','a');
        fidCplex = fopen('logCplex.txt','a');
        
    end
    fprintf(fid,'********************* Iteration %3.d ********************\n',opt.it);
    fprintf(fid,'it xi      slack   AddCons LP-TIME ST infViol NewCons QP-TIME ST\n');

    
    fprintf(fidCplex,'********************* Iteration %3.d ********************\n',opt.it);
    
    DisplayFunc = @(x) fprintf(fidCplex,[x,'\n']);
else
    DisplayFunc = [];
end



uDims = cellfun(@numel,ldu);
xDims = cellfun(@(z)cellfun(@numel,z),ldx,'UniformOutput',false);
vDims = cellfun(@(z)cellfun(@numel,z),ldv,'UniformOutput',false);
sDims = numel(lds);


nuH = sum(uDims);


du = [];
dx = [];
dv = [];
ds = [];

violationH = [];


lowActiveP = opt.lowActive;
upActiveP = opt.upActive;

linesAdded.ub.x = upActiveP.x;
linesAdded.lb.x = lowActiveP.x;
linesAdded.ub.v = upActiveP.v;
linesAdded.lb.v = lowActiveP.v;
linesAdded.ub.s = upActiveP.s;
linesAdded.lb.s = lowActiveP.s;


newCons = cell(opt.maxQpIt+1,1);
nc = cell(opt.maxQpIt+1,1);

% mantain a record of the constraints added at each iteration
newCons{1} = linesAdded;

% count the number of new constraints
nc{1}.ub.x = cellfun(@(act)sum(cellfun(@sum,act)),linesAdded.ub.x,'UniformOutput',false);
nc{1}.lb.x = cellfun(@(act)sum(cellfun(@sum,act)),linesAdded.lb.x,'UniformOutput',false);
nc{1}.ub.v = cellfun(@(act)sum(cellfun(@sum,act)),linesAdded.ub.v,'UniformOutput',false);
nc{1}.lb.v = cellfun(@(act)sum(cellfun(@sum,act)),linesAdded.lb.v,'UniformOutput',false);
nc{1}.ub.s = cellfun(@sum,linesAdded.ub.s);
nc{1}.lb.s = cellfun(@sum,linesAdded.lb.s);

violation = [];
solved = false;


nAddRows = 0;


P = Cplex('LP - QP');
P.DisplayFunc = DisplayFunc;
P.Param.qpmethod.Cur = 6;
P.Param.lpmethod.Cur = 6;
P.Param.emphasis.numerical.Cur = 1;


% Add control variables
P.addCols(zeros(nuH+2,1), [], [0;0;cell2mat(ldu)], [1;1/opt.bigM;cell2mat(udu)]);  % objective will be set up again later
P.Param.timelimit.Cur = 7200;


for k = 1:opt.maxQpIt
    
    
    if opt.condense
    
        % extract the current constraints lines to add in this iteration
        [ Alow.x,blow.x,xRlow ] = cellfun(@(A,ld,newConsrk,a)buildActiveConstraints(A,ld,newConsrk,-1,'R',a),Ax,ldx,newCons{k}.lb.x,ax,'UniformOutput',false);
        [ Aup.x,bup.x,xRup ] =    cellfun(@(A,ld,newConsrk,a)buildActiveConstraints(A,ld,newConsrk, 1,'R',a),Ax,udx,newCons{k}.ub.x,ax,'UniformOutput',false);
        [ Alow.v,blow.v,vRlow ] = cellfun(@(A,ld,newConsrk,a)buildActiveConstraints(A,ld,newConsrk,-1,'R',a),Av,ldv,newCons{k}.lb.v,av,'UniformOutput',false);
        [ Aup.v,bup.v,vRup ] =    cellfun(@(A,ld,newConsrk,a)buildActiveConstraints(A,ld,newConsrk, 1,'R',a),Av,udv,newCons{k}.ub.v,av,'UniformOutput',false);       
        
        Aup.s = cellfun(@(Asu)Asu(cell2mat(newCons{k}.ub.s),:),As,'UniformOutput',false);
        bup.s = uds(cell2mat(newCons{k}.ub.s));
        sRup = as(cell2mat(newCons{k}.ub.s));        
        
        Alow.s = cellfun(@(Asu)-Asu(cell2mat(newCons{k}.lb.s),:),As,'UniformOutput',false);
        blow.s = -uds(cell2mat(newCons{k}.lb.s));
        sRlow = -as(cell2mat(newCons{k}.lb.s));
            
        Aact = cell2mat([cellfun(@(A)cell2mat(A),Alow.x,'UniformOutput',false);
                         cellfun(@(A)cell2mat(A), Aup.x,'UniformOutput',false);
                         cellfun(@(A)cell2mat(A),Alow.v,'UniformOutput',false);
                         cellfun(@(A)cell2mat(A), Aup.v,'UniformOutput',false);
                         {cell2mat(Alow.s)};
                         {cell2mat(Aup.s)}]);
              
    else
        blow.x = cellfun(@(xr,ar)cellfun(@(xrk,ark)-xrk(ark),xr ,ar,'UniformOutput',false),ldx,newCons{k}.lb.x,'UniformOutput',false);
        xRlow =  cellfun(@(xr,ar)cellfun(@(xrk,ark)-xrk(ark),xr ,ar,'UniformOutput',false),ax ,newCons{k}.lb.x,'UniformOutput',false);
        bup.x =  cellfun(@(xr,ar)cellfun(@(xrk,ark) xrk(ark),xr ,ar,'UniformOutput',false),udx,newCons{k}.ub.x,'UniformOutput',false);
        xRup =   cellfun(@(xr,ar)cellfun(@(xrk,ark) xrk(ark),xr ,ar,'UniformOutput',false),ax ,newCons{k}.ub.x,'UniformOutput',false);
        blow.v = cellfun(@(xr,ar)cellfun(@(xrk,ark)-xrk(ark),xr ,ar,'UniformOutput',false),ldv,newCons{k}.lb.v,'UniformOutput',false);
        vRlow =  cellfun(@(xr,ar)cellfun(@(xrk,ark)-xrk(ark),xr ,ar,'UniformOutput',false),av ,newCons{k}.lb.v,'UniformOutput',false);
        bup.v =  cellfun(@(xr,ar)cellfun(@(xrk,ark) xrk(ark),xr ,ar,'UniformOutput',false),udv,newCons{k}.ub.v,'UniformOutput',false);
        vRup =   cellfun(@(xr,ar)cellfun(@(xrk,ark) xrk(ark),xr ,ar,'UniformOutput',false),av ,newCons{k}.ub.v,'UniformOutput',false);

        bup.s = uds(cell2mat(newCons{k}.ub.s));
        sRup = as(cell2mat(newCons{k}.ub.s));        
        
        blow.s = -uds(cell2mat(newCons{k}.lb.s));
        sRlow = -as(cell2mat(newCons{k}.lb.s));        
        
        
        if (k == 1 && ~isempty(Aact1))
            Aact = cell2mat(Aact1);
        else
            Aact = constraintBuilder(newCons{k});
            Aact = cell2mat(Aact{1});
        end
        
    end
    

    xRs = cell2mat([[
        cellfun(@(R)cell2mat(R),xRlow,'UniformOutput',false);
        cellfun(@(R)cell2mat(R),xRup,'UniformOutput',false);
        cellfun(@(R)cell2mat(R),vRlow,'UniformOutput',false);
        cellfun(@(R)cell2mat(R),vRup,'UniformOutput',false);
        {sRlow};
        {sRup }], ...
        [cellfun(@(xir)cell2mat(cellfun(@(xirk)-ones(sum(xirk),1),xir,'UniformOutput',false)),newCons{k}.lb.x,'UniformOutput',false);...
         cellfun(@(xir)cell2mat(cellfun(@(xirk)-ones(sum(xirk),1),xir,'UniformOutput',false)),newCons{k}.ub.x,'UniformOutput',false);...
         cellfun(@(xir)cell2mat(cellfun(@(xirk)-ones(sum(xirk),1),xir,'UniformOutput',false)),newCons{k}.lb.v,'UniformOutput',false);...
         cellfun(@(xir)cell2mat(cellfun(@(xirk)-ones(sum(xirk),1),xir,'UniformOutput',false)),newCons{k}.ub.v,'UniformOutput',false);...
         {-ones(sum(cell2mat(newCons{k}.lb.s)),1)};
         {-ones(sum(cell2mat(newCons{k}.ub.s)),1)}]]);           

 
                
    % Merge constraints in one matrix !
    
    Aq = [xRs,Aact];
    bq = cell2mat([
    cellfun(@(R)cell2mat(R),blow.x,'UniformOutput',false);
	cellfun(@(R)cell2mat(R),bup.x,'UniformOutput',false);    
    cellfun(@(R)cell2mat(R),blow.v,'UniformOutput',false);    
    cellfun(@(R)cell2mat(R),bup.v,'UniformOutput',false);    
    {blow.s}
    {bup.s}]);
   
    
    if opt.qpDebug
        fprintf(fidCplex,'************* LP %d  ******************\n',k);
    end
    
    
    % now add the new rows
    P.addRows(-inf(size(bq)),Aq,bq);
    
	P.Model.lb(1) = 0;
    P.Model.ub(1) = 1;
	P.Model.lb(2) = 0;
    P.Model.ub(2) = 1/opt.bigM;
    
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
        fprintf(fid,'%2.d %1.1e %1.1e %1.1e %1.1e %2.d ',k,1-xibar,sM,nAddRows,lpTime+lpTime2,P.Solution.status) ;
        fprintf(fidCplex,'************* QP %d  *******************\n',k);
    end
    
    % set up the qp objective
    P.Model.Q = blkdiag(1,1,M);
    B =cell2mat(g) + xibar * cell2mat(w);
    P.Model.obj = [0;0;B'];
    P.Model.lb(1) = xibar;
    P.Model.ub(1) = xibar;
    P.Model.lb(2) = sM;
    P.Model.ub(2) = sM;
    
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
    
    
    duC = mat2cell(du,uDims,1);
    
    if opt.condense

        dxN = cellfun(@(Ar,ssr)cellmtimes(Ar,duC,'lowerTriangular',true,'ci',ssr.ci),Ax,opt.ss,'UniformOutput',false);
        dvN = cellfun(@(Ar,ssr)cellmtimes(Ar,duC,'lowerTriangular',true,'ci',ssr.ci),Av,opt.ss,'UniformOutput',false);
        dsN = cell2mat(As)*cell2mat(duC);

    
    else
        
        [dxN,dvN,dsN] = predictor(duC);
    end
    
    du = duC;
    dx = cellfun(@(ar,dr)cellfun(@(z,dz)xibar*z+dz,ar,dr,'UniformOutput',false),ax,dxN,'UniformOutput',false);
	dv = cellfun(@(ar,dr)cellfun(@(z,dz)xibar*z+dz,ar,dr,'UniformOutput',false),av,dvN,'UniformOutput',false);
    ds = as*xibar + dsN;

    
    % Check which other constraints are infeasible
    [feasible.x,lowActive.x,upActive.x,violation.x ] = cellfun(@(d,l,u)checkConstraintFeasibility(d,l,u,'primalFeasTol',opt.feasTol ),dx,ldx,udx,'UniformOutput',false) ;
    [feasible.v,lowActive.v,upActive.v,violation.v ] = cellfun(@(d,l,u)checkConstraintFeasibility(d,l,u,'primalFeasTol',opt.feasTol ),dv,ldv,udv,'UniformOutput',false) ;
    [feasible.s,lowActive.s,upActive.s,violation.s ] = checkConstraintFeasibility({ds},{lds},{uds},'primalFeasTol',opt.feasTol );

    violation.x = max(cell2mat(violation.x));
    violation.v = max(cell2mat(violation.v));
    
    % debugging purpouse:  see if the violation is decreasing!
    ineqViolation = violation.x;
    ineqViolation = max(ineqViolation,violation.v);
    ineqViolation = max(ineqViolation,violation.s);
    
    
    violationH = [violationH,ineqViolation];
    
    
    
    % Determine the new set of contraints to be added
    [newCons{k+1}.lb.x,nc{k+1}.lb.x] =  cellfun(@(b,n)newConstraints(b,n),linesAdded.lb.x, lowActive.x,'UniformOutput',false);
    [newCons{k+1}.ub.x,nc{k+1}.ub.x] =  cellfun(@(b,n)newConstraints(b,n),linesAdded.ub.x, upActive.x,'UniformOutput',false);
	[newCons{k+1}.lb.v,nc{k+1}.lb.v] =  cellfun(@(b,n)newConstraints(b,n),linesAdded.lb.v, lowActive.v,'UniformOutput',false);
    [newCons{k+1}.ub.v,nc{k+1}.ub.v] =  cellfun(@(b,n)newConstraints(b,n),linesAdded.ub.v, upActive.v,'UniformOutput',false);
    [newCons{k+1}.lb.s,nc{k+1}.lb.s] =  newConstraints(linesAdded.lb.s, lowActive.s);
    [newCons{k+1}.ub.s,nc{k+1}.ub.s] =  newConstraints(linesAdded.ub.s, upActive.s);
    
    % determine how many new constraints are being added in total
    newC = 0;
    bN = {'lb','ub'};
    vN = {'x','v'};
    for bi = 1:numel(bN)
        for vi = 1:numel(vN);
            newC  = newC + sum(cell2mat(nc{k+1}.(bN{bi}).(vN{vi})));
        end
    end
    newC = newC + nc{k+1}.ub.s + nc{k+1}.lb.s;
    
    if opt.qpDebug
        fprintf(fid,'%1.1e %1.1e %1.1e %2.d\n',ineqViolation,newC,QpTime+QpTime2,P.Solution.status) ;
    end
    
    % if we cannot add more constraints, so the problem is solved!
    if newC == 0
        if ineqViolation > opt.feasTol
            if opt.qpDebug
                fprintf(fid,'Irreductible constraint violation inf norm: %e \n',ineqViolation) ;
            end
        end
        solved = true;
        break;
    end
    
    
    linesAdded.lb.x = cellfun(@(a1,a2)unionActiveSets(a1,a2),linesAdded.lb.x,lowActive.x,'UniformOutput',false);
    linesAdded.ub.x = cellfun(@(a1,a2)unionActiveSets(a1,a2),linesAdded.ub.x,upActive.x,'UniformOutput',false);
	linesAdded.lb.v = cellfun(@(a1,a2)unionActiveSets(a1,a2),linesAdded.lb.v,lowActive.v,'UniformOutput',false);
	linesAdded.ub.v = cellfun(@(a1,a2)unionActiveSets(a1,a2),linesAdded.ub.v,upActive.v,'UniformOutput',false);
	linesAdded.lb.s = unionActiveSets(linesAdded.lb.s,lowActive.s);
	linesAdded.ub.s = unionActiveSets(linesAdded.ub.s,upActive.s);
    
end
if ~solved
    warning('Qps were not solved within the iteration limit')
end

xi = 1-xibar;
slack = sM;

if isfield(P.Solution,'dual')
    dual = P.Solution.dual;
else
    dual = [];
end


% extract the dual variables:
to = 0;  %% skip xibar and sM dual constraint!
mu.lb.x = cell(1,numel(xDims));
mu.ub.x = cell(1,numel(xDims));
mu.lb.v = cell(1,numel(vDims));
mu.ub.v = cell(1,numel(vDims));

returnRows = true;
for j = 1:k
    
    for r = 1:numel(xDims)
        from = to + 1;
        to = from + nc{j}.lb.x{r}-1;
        [d] = extractCompressIneq(-dual(from:to),newCons{j}.lb.x{r},xDims{r},returnRows);
        if j==1
            mu.lb.x{r} = d;
        else
            mu.lb.x{r} = cellfun(@(x1,x2)x1+x2,mu.lb.x{r},d,'UniformOutput',false);
        end
    end
    
    for r = 1:numel(xDims)
        from = to + 1;
        to = from + nc{j}.ub.x{r}-1;
        [d] = extractCompressIneq(-dual(from:to),newCons{j}.ub.x{r},xDims{r},returnRows);
        if j==1
            mu.ub.x{r} = d;
        else
            mu.ub.x{r} = cellfun(@(x1,x2)x1+x2,mu.ub.x{r},d,'UniformOutput',false);
        end
    end
    
    for r = 1:numel(vDims)
        from = to + 1;
        to = from + nc{j}.lb.v{r}-1;
        [d] = extractCompressIneq(-dual(from:to),newCons{j}.lb.v{r},vDims{r},returnRows);
        if j==1
            mu.lb.v{r} = d;
        else
            mu.lb.v{r} = cellfun(@(x1,x2)x1+x2,mu.lb.v{r},d,'UniformOutput',false);
        end
    end
    
    for r = 1:numel(vDims)
        from = to + 1;
        to = from + nc{j}.ub.v{r}-1;
        [d] = extractCompressIneq(-dual(from:to),newCons{j}.ub.v{r},vDims{r},returnRows);
        if j==1
            mu.ub.v{r} = d;
        else
            mu.ub.v{r} = cellfun(@(x1,x2)x1+x2,mu.ub.v{r},d,'UniformOutput',false);
        end
    end   

    from = to + 1;
    to = from + nc{j}.lb.s-1;
    [d] = extractCompressIneq(-dual(from:to),newCons{j}.lb.s,sDims,returnRows);
    if j==1
        mu.lb.s = d;
    else
        mu.lb.s = cellfun(@(x1,x2)x1+x2,mu.lb.s,d,'UniformOutput',false);
    end


    from = to + 1;
    to = from + nc{j}.ub.s-1;
    [d] = extractCompressIneq(-dual(from:to),newCons{j}.ub.s,sDims,returnRows);
    if j==1
        mu.ub.s = d;
    else
        mu.ub.s = cellfun(@(x1,x2)x1+x2,mu.ub.s,d,'UniformOutput',false);
    end

    
end
if to ~= size(P.Model.rhs,1)
    error('checkSizes!');
end

%%%  First order optimality
%norm(P.Model.Q * P.Solution.x + P.Model.obj -(P.Model.A)' * dual - P.Solution.reducedcost)


% extract dual variables with respect to the controls
mu.ub.u = mat2cell(max(-P.Solution.reducedcost(3:nuH+2),0)',1,uDims);
mu.lb.u = mat2cell(max( P.Solution.reducedcost(3:nuH+2),0)',1,uDims);

if opt.qpDebug
    if opt.condense
        optCheck = cell2mat(duC)'*M + B + ...
            catAndSum(cellfun(@(mup,mul,A,ssr)cell2mat(cellmtimesT( cellfun(@(x1,x2)x1-x2,mup,mul,'UniformOutput',false),A ,'lowerTriangular',true,'ci',ssr.ci,'columnVector',false)),mu.ub.x',mu.lb.x',Ax,opt.ss,'UniformOutput',false)) + ...
            catAndSum(cellfun(@(mup,mul,A,ssr)cell2mat(cellmtimesT( cellfun(@(x1,x2)x1-x2,mup,mul,'UniformOutput',false),A ,'lowerTriangular',true,'ci',ssr.ci,'columnVector',false)),mu.ub.v',mu.lb.v',Av,opt.ss,'UniformOutput',false)) + ...
            (cell2mat(mu.ub.s)-cell2mat(mu.lb.s))*cell2mat(As) +...
            cell2mat(             cellfun(@(x1,x2)x1-x2,mu.ub.u,mu.lb.u,'UniformOutput',false));


        optNorm = norm(optCheck);
        fprintf(fid,'Optimality norm: %e \n',optNorm) ;
        if optNorm > opt.feasTol*10
            warning('QP optimality norm might be to high');
        end
        
    elseif opt.testQP
        
        J.Jx = minusC(mu.ub.x,mu.lb.x);
        J.Jv = minusC(mu.ub.v,mu.lb.v);
        J.Ju = minusS(mu.ub.u,mu.lb.u);
        J.Js = minusS(mu.ub.s,mu.lb.s);
        J.Js = cell2mat(J.Js);
                
        optCheck = cell2mat(duC)'*M  + B + cell2mat(opt.lagFunc(J)) ;
        optNorm = norm(optCheck);
        fprintf(fid,'Optimality norm: %e \n',optNorm) ;
        if optNorm > opt.feasTol*10
            warning('QP optimality norm might be to high');
        end
        
    end
end


qpVAl = P.Solution.objval;

if opt.qpDebug
    
    
    fclose(fid);
end

end


function out = catAndSum(M)


if any(cellfun(@issparse,M))
    if isrow(M)
        M = M';
    end
    rows= size(M{1},1);
    blocks = numel(M);
    out = sparse( repmat(1:rows,1,blocks),1:rows*blocks,1)*cell2mat(M);
else
    out = sum(cat(3,M{:}),3);    
end


end


function me = minusC(mU,mL)
    f = @minusS;
    me = cellfun(f,mU,mL,'UniformOutput',false);
end

function me = minusS(mU,mL)
    me = cellfun(@(x1,x2)(x1-x2),mU,mL,'UniformOutput',false);
end