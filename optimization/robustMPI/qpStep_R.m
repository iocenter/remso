function [ duC,dx,dv,ds,xi,lowActive,upActive,mu,violation,qpVAl,dxN,dvN,dsN,slack,k] = qpStep_R(M,g,w,ldu,udu,Aact1,predictor,constraintBuilder,ax,Ax,ldx,udx,av,Av,ldv,udv,as,As,lds,uds,sss,varargin )
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


opt = struct('qpDebug',true,'lowActive',[],'upActive',[],'feasTol',1e-6,'maxQpIt',20,'it',0,'bigM',1e9,'condense',true,'lagFunc',[],'testQP',false,'nCons',100);
opt = merge_options(opt, varargin{:});

ss = sss.ss;
jobSchedule = sss.jobSchedule;
imMaster = jobSchedule.imMaster;

%assert(opt.feasTol*10>1/opt.bigM,'Make sure bigM is big enough compared to feasTol');
opt.bigM = max(opt.bigM,10/opt.feasTol);
xibar = 1;
feasTol = opt.feasTol;
nCons = opt.nCons;

if opt.qpDebug && imMaster
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
%spmd
xDims = getZDims(ldx);
vDims = getZDims(ldv);
%end
sDims = numel(lds);


nuH = sum(uDims);


du = [];
dx = [];
dv = [];
ds = [];

violationH = [];


lowActiveP = opt.lowActive;
upActiveP = opt.upActive;

linesAddedlbx = lowActiveP.x;
linesAddedubx = upActiveP.x;
linesAddedlbv = lowActiveP.v;
linesAddedubv = upActiveP.v;
if imMaster
linesAddedlbs = lowActiveP.s;
linesAddedubs = upActiveP.s;
end

newlbx = linesAddedlbx;
newubx = linesAddedubx;
newlbv = linesAddedlbv;
newubv = linesAddedubv;
if imMaster
newlbs = linesAddedlbs;
newubs = linesAddedubs;
end


newCons = cell(opt.maxQpIt+1,1);
nc = cell(opt.maxQpIt+1,1);

% mantain a record of the constraints added at each iteration
newCons{1}.lb.x = newlbx;
newCons{1}.ub.x = newubx;
newCons{1}.lb.v = newlbv;
newCons{1}.ub.v = newubv;
if imMaster
newCons{1}.lb.s = newlbs;
newCons{1}.ub.s = newubs;
end

% count the number of new constraints
%spmd
nclbx = countActiveLines(newlbx);
ncubx = countActiveLines(newubx);
nclbv = countActiveLines(newlbv);
ncubv = countActiveLines(newubv);
%end
if imMaster
nclbs = cellfun(@sum,newlbs);
ncubs = cellfun(@sum,newubs);
end


% count the number of new constraints
nc{1}.lb.x = nclbx;
nc{1}.ub.x = ncubx;
nc{1}.lb.v = nclbv;
nc{1}.ub.v = ncubv;
if imMaster
nc{1}.lb.s = nclbs;
nc{1}.ub.s = ncubs;
end

violation = [];
solved = false;


nAddRows = 0;

if imMaster
P = Cplex('LP - QP');
P.DisplayFunc = DisplayFunc;
P.Param.qpmethod.Cur = 6;
P.Param.lpmethod.Cur = 6;
P.Param.emphasis.numerical.Cur = 1;


% Add control variables
P.addCols(zeros(nuH+2,1), [], [0;0;cell2mat(ldu)], [1;1/opt.bigM;cell2mat(udu)]);  % objective will be set up again later
P.Param.timelimit.Cur = 7200;

Q = blkdiag(1,1,M);
% CPLEX keeps complaining that the approximation is not symmetric
if ~issymmetric(Q)
	warning('Hessian approximation seems to be not symmetric')
end
Q = (Q+Q')/2;
end

for k = 1:opt.maxQpIt
    
    
    if opt.condense
    
        % extract the current constraints lines to add in this iteration
        %spmd
        [Alowx,blowx,xRlow] = applyBuildConstraints(-1,Ax,ldx,newlbx,ax);
        [Aupx ,bupx ,xRup ] = applyBuildConstraints( 1,Ax,udx,newubx,ax);
        [Alowv,blowv,vRlow] = applyBuildConstraints(-1,Av,ldv,newlbv,av);
        [Aupv ,bupv ,vRup ] = applyBuildConstraints( 1,Av,udv,newubv,av);
        %end  
        
        error('not Implemented: bringVariablesMPI should cend a matrix!')
        Alowx = bringVariablesMPI(Alowx,jobSchedule);
        blowx = bringVariablesMPI(blowx,jobSchedule);
        xRlow = bringVariablesMPI(xRlow,jobSchedule);
        
        Aupx =  bringVariablesMPI(Aupx,jobSchedule);
        bupx =  bringVariablesMPI(bupx,jobSchedule);
        xRup =  bringVariablesMPI(xRup,jobSchedule);
        
        Alowv = bringVariablesMPI(Alowv,jobSchedule);
        blowv = bringVariablesMPI(blowv,jobSchedule);
        vRlow = bringVariablesMPI(vRlow,jobSchedule);
        
        Aupv =  bringVariablesMPI(Aupv,jobSchedule);
        bupv =  bringVariablesMPI(bupv,jobSchedule);
        vRup =  bringVariablesMPI(vRup,jobSchedule);
        
               
        if imMaster
        Aups = cellfun(@(Asu)Asu(cell2mat(newubs),:),As,'UniformOutput',false);
        bups = uds(cell2mat(newubs));
        sRup = as(cell2mat(newubs));        
        
        Alows = cellfun(@(Asu)-Asu(cell2mat(newlbs),:),As,'UniformOutput',false);
        blows = -uds(cell2mat(newlbs));
        sRlow = -as(cell2mat(newlbs));
            
        Aact = cell2mat(vertcat(...
        vertcat(Alowx{:}),...
        vertcat(Aupx{:}),...
        vertcat(Alowv{:}),...
        vertcat(Aupv{:}),...
        vertcat(Alows{:}),...
        vertcat(Aups{:})));
        end
        
              
    else
        %spmd
        blowx = extractActiveVector(ldx,newlbx,-1);
        xRlow = extractActiveVector(ax ,newlbx,-1);
        bupx  = extractActiveVector(udx,newubx, 1);
        xRup  = extractActiveVector(ax ,newubx, 1);
        blowv = extractActiveVector(ldv,newlbv,-1);
        vRlow = extractActiveVector(av ,newlbv,-1);
        bupv  = extractActiveVector(udv,newubv, 1);
        vRup  = extractActiveVector(av ,newubv, 1);
        %end
   
        blowx = bringVariablesMPI(blowx,jobSchedule);
        xRlow = bringVariablesMPI(xRlow,jobSchedule);
        
        bupx =  bringVariablesMPI(bupx,jobSchedule);
        xRup =  bringVariablesMPI(xRup,jobSchedule);
        
        blowv = bringVariablesMPI(blowv,jobSchedule);
        vRlow = bringVariablesMPI(vRlow,jobSchedule);
        
        bupv =  bringVariablesMPI(bupv,jobSchedule);
        vRup =  bringVariablesMPI(vRup,jobSchedule);        
        
        
        if imMaster
        bups = uds(cell2mat(newubs));
        sRup = as(cell2mat(newubs));        
        
        blows = -uds(cell2mat(newlbs));
        sRlow = -as(cell2mat(newlbs));        
        end
       
        if (k == 1 && ~isempty(Aact1))
            Aact = cell2mat(Aact1);
        else
            Aact = constraintBuilder(newCons{k});
            Aact = cell2mat(Aact{1});
        end
        
    end
    
    if imMaster 
    xRs = cell2mat(vertcat(...
        vertcat(xRlow{:}),...
        vertcat(xRup{:}),...
        vertcat(vRlow{:}),...
        vertcat(vRup{:}),...
        {sRlow},...
        {sRup }));
    xRs = [xRs,-ones(size(xRs,1),1)];
        
 
                
    % Merge constraints in one matrix !
    
    Aq = [xRs,Aact];
    bq = cell2mat(vertcat(...
    blowx{:},...
	bupx{:},...   
    blowv{:},...    
    bupv{:},...    
    {blows},...
    {bups}));
   
    
    if opt.qpDebug && imMaster
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
            if opt.qpDebug && imMaster
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
    xibar = min(max(P.Solution.x(1),0),1);
    sM = max(P.Solution.x(2),0);
    
    if opt.qpDebug && imMaster
        fprintf(fid,'%2.d %1.1e %1.1e %1.1e %1.1e %2.d ',k,1-xibar,sM,nAddRows,lpTime+lpTime2,P.Solution.status) ;
        fprintf(fidCplex,'************* QP %d  *******************\n',k);
    end
    
    % set up the qp objective
    err = checkSolutionFeasibility(P);
    P.Model.Q = Q;
    B =cell2mat(g) + xibar * cell2mat(w);
    P.Model.obj = [0;opt.bigM;B'];
    P.Model.lb(1) = xibar;
    P.Model.ub(1) = xibar;
    %P.Model.lb(2) = sM;
    P.Model.ub(2) = sM + err;
    
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
            if opt.qpDebug && imMaster
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
    else
	du = zeros(nuH,1);
    end
    du=NMPI_Bcast(du,length(du),jobSchedule.Master_rank,jobSchedule.my_rank); 
    
    duC = mat2cell(du,uDims,1);
    
    if opt.condense
        %spmd
        dxN = nullSpaceStep(Ax,duC,ss);
        dvN = nullSpaceStep(Av,duC,ss);
        %end
		if imMaster
        dsN = cell2mat(As)*cell2mat(duC);
		end
    else
        [dxN,dvN,dsN] = predictor(duC);
    end
    
    du = duC;
    %spmd
    dx = composeStep(ax,dxN,xibar);
	dv = composeStep(av,dvN,xibar);
    %end

    if imMaster
    ds = as*xibar + dsN;
    end

    %spmd
    % Check which other constraints are infeasible
    [feasiblex,lowActivex,upActivex,violationx ] = applyCheckConstraintFeasibility(dx,ldx,udx,0,nCons)  ;
    [feasiblev,lowActivev,upActivev,violationv ] = applyCheckConstraintFeasibility(dv,ldv,udv,0,nCons)  ;
    violationx = max([violationx;-inf]);
    violationv = max([violationv;-inf]);
    violationx = gopMPI('M',violationx,jobSchedule);
    violationv = gopMPI('M',violationv,jobSchedule);
	violationxv = NMPI_Bcast([violationx,violationv],2,jobSchedule.Master_rank,jobSchedule.my_rank);  %% not extrictly necessary but done to mantain same results with the othe algs
    %end
    
    if imMaster
    [feasibles,lowActives,upActives,violations ] = checkConstraintFeasibility({ds},{lds},{uds},'primalFeasTol',0,'first',nCons);
    end
    violationx = violationxv(1);
    violationv = violationxv(2);
    
    % debugging purpouse:  see if the violation is decreasing!
    ineqViolation = violationx;
    ineqViolation = max(ineqViolation,violationv);
    if imMaster
    ineqViolation = max(ineqViolation,violations);
    
    
    violationH = [violationH,ineqViolation];
    end
    
    
    % Determine the new set of contraints to be added
    %spmd
    [newlbx,nclbx] =  cellfun(@newConstraints,linesAddedlbx, lowActivex,'UniformOutput',false);
    [newubx,ncubx] =  cellfun(@newConstraints,linesAddedubx, upActivex,'UniformOutput',false);
	[newlbv,nclbv] =  cellfun(@newConstraints,linesAddedlbv, lowActivev,'UniformOutput',false);
    [newubv,ncubv] =  cellfun(@newConstraints,linesAddedubv, upActivev,'UniformOutput',false);
    nclbxCount = sum([cell2mat(nclbx);0]);
    ncubxCount = sum([cell2mat(ncubx);0]);
    nclbvCount = sum([cell2mat(nclbv);0]);
    ncubvCount = sum([cell2mat(ncubv);0]);
    
    nclbxCount = gopMPI('+',nclbxCount,jobSchedule);
    ncubxCount = gopMPI('+',ncubxCount,jobSchedule);
    nclbvCount = gopMPI('+',nclbvCount,jobSchedule);
    ncubvCount = gopMPI('+',ncubvCount,jobSchedule);
    
    %end
    if imMaster
    [newlbs,nclbs] =  newConstraints(linesAddedlbs, lowActives);
    [newubs,ncubs] =  newConstraints(linesAddedubs, upActives);
	nclbsCount = nclbs;
    ncubsCount = ncubs;
    
    newC = nclbxCount+ncubxCount+nclbvCount+ncubvCount+nclbsCount+ncubsCount;
    end
    
    newCons{k+1}.lb.x = newlbx;
    newCons{k+1}.ub.x = newubx;
    newCons{k+1}.lb.v = newlbv;
    newCons{k+1}.ub.v = newubv;
	if imMaster
    newCons{k+1}.lb.s = newlbs;
    newCons{k+1}.ub.s = newubs;
	end    
    nc{k+1}.lb.x = nclbx;
    nc{k+1}.ub.x = ncubx;
    nc{k+1}.lb.v = nclbv;
    nc{k+1}.ub.v = ncubv;
	if imMaster
    nc{k+1}.lb.s = nclbs;
    nc{k+1}.ub.s = ncubs;
    
    if opt.qpDebug && imMaster
        fprintf(fid,'%1.1e %1.1e %1.1e %2.d\n',ineqViolation,newC,QpTime+QpTime2,P.Solution.status) ;
    end
    
    % if we cannot add more constraints, so the problem is solved!
    solved = false;
    if newC == 0 || ineqViolation < feasTol
        if ineqViolation > feasTol
            if opt.qpDebug && imMaster
                fprintf(fid,'Irreductible constraint violation inf norm: %e \n',ineqViolation) ;
            end
        end
        solved = true;
    end
	end
	solved=NMPI_Bcast(solved+0,length(solved),jobSchedule.Master_rank,jobSchedule.my_rank)==1;
	if solved
		break;
	end
    %spmd
    linesAddedlbx = cellfun(@unionActiveSets,linesAddedlbx,lowActivex,'UniformOutput',false);
    linesAddedubx = cellfun(@unionActiveSets,linesAddedubx,upActivex,'UniformOutput',false);
	linesAddedlbv = cellfun(@unionActiveSets,linesAddedlbv,lowActivev,'UniformOutput',false);
	linesAddedubv = cellfun(@unionActiveSets,linesAddedubv,upActivev,'UniformOutput',false);
    %end
	if imMaster
	linesAddedlbs = unionActiveSets(linesAddedlbs,lowActives);
	linesAddedubs = unionActiveSets(linesAddedubs,upActives);
	end
end
if ~solved
    warning('Qps were not solved within the iteration limit')
end
if imMaster
xi = 1-xibar;
slack = sM;

if isfield(P.Solution,'dual')
    dual = P.Solution.dual;
else
    dual = [];
end

else
dual =[];
end
% extract the dual variables:
to = 0;  %% skip xibar and sM dual constraint!
mu.dx = [];
mu.dv = [];
mu.ds = [];


for j = 1:k
    
    first = j==1;
    
    [mu.dx,to] = extractDuals(dual,to,nc{j}.lb.x,newCons{j}.lb.x,xDims,mu.dx,first,-1,jobSchedule);    
	[mu.dx,to] = extractDuals(dual,to,nc{j}.ub.x,newCons{j}.ub.x,xDims,mu.dx,false, 1,jobSchedule);
    
    [mu.dv,to] = extractDuals(dual,to,nc{j}.lb.v,newCons{j}.lb.v,vDims,mu.dv,first,-1,jobSchedule);   
	[mu.dv,to] = extractDuals(dual,to,nc{j}.ub.v,newCons{j}.ub.v,vDims,mu.dv,false, 1,jobSchedule);
    
	if imMaster
    [mu.ds,to] = extractDualsS(dual,to,nc{j}.lb.s,newCons{j}.lb.s,sDims,mu.ds,first,-1);    
	[mu.ds,to] = extractDualsS(dual,to,nc{j}.ub.s,newCons{j}.ub.s,sDims,mu.ds,false, 1);   
	end


    
end
if imMaster && to ~= size(P.Model.rhs,1) 
    error('checkSizes!');
end

%%%  First order optimality
%norm(P.Model.Q * P.Solution.x + P.Model.obj -(P.Model.A)' * dual - P.Solution.reducedcost)


% extract dual variables with respect to the controls
if imMaster
mudu = P.Solution.reducedcost(3:nuH+2);
else
mudu = zeros(nuH,1);
end
mudu = NMPI_Bcast(mudu,length(mudu),jobSchedule.Master_rank,jobSchedule.my_rank);
mu.du = mat2cell(-mudu',1,uDims);

if opt.qpDebug
    if opt.condense

        mudx = mu.dx;
        mudv = mu.dv;
        %spmd
           optCheck = gradientLagrangian(mudx,Ax,ss) + gradientLagrangian(mudv,Av,ss);
           optCheck = gopMPI('+',optCheck,jobSchedule);
        %end

        if imMaster
        optCheck = cell2mat(duC)'*M + B + optCheck + (cell2mat(mu.ds))*cell2mat(As) + cell2mat(mu.du);
        
        optNorm = norm(optCheck);
        fprintf(fid,'Optimality norm: %e \n',optNorm) ;
        if optNorm > opt.feasTol*10
            warning('QP optimality norm might be to high');
        end
        end
    elseif opt.testQP
        
        J.Jx = mu.dx;
        J.Jv = mu.dv;
        J.Ju = mu.du;
        if imMaster
        J.Js = mu.ds;
        J.Js = cell2mat(J.Js);
        else
        J.Js = nan;
        end        
        optCheck = cell2mat(opt.lagFunc(J)) ;
        if imMaster
        optCheck = cell2mat(duC)'*M  + B + optCheck;
        optNorm = norm(optCheck);
        fprintf(fid,'Optimality norm: %e \n',optNorm) ;
        if optNorm > opt.feasTol*10
            warning('QP optimality norm might be to high');
        end
        end
    end
end

% Make sure that only the non-weakly active constraints are kept for the
% next iteration.
mudx = mu.dx;
mudv = mu.dv;
%spmd
lowActivex = getLowerActives(mudx);
upActivex =  getUpperActives(mudx);
lowActivev = getLowerActives(mudv);
upActivev =  getUpperActives(mudv);
%end
lowActive.x = lowActivex;
upActive.x = upActivex;
lowActive.v =lowActivev;
upActive.v = upActivev;

if imMaster
lowActive.s = cellfun(@(l)l'<0,mu.ds,'UniformOutput',false);
upActive.s =  cellfun(@(l)l'>0,mu.ds,'UniformOutput',false);
end
violation.x = violationx;
violation.v = violationv;
if imMaster
violation.s = violations;


qpVAl = P.Solution.objval;


if opt.qpDebug
    
    
    fclose(fid);
end

else
qpVAl = 0;
xi = -1;
slack = -1;
dsN = [];
end
xi = NMPI_Bcast(xi,1,jobSchedule.Master_rank,jobSchedule.my_rank);


end

function [a] = getLowerActives(mudz)
a = cellfun(@(lr)cellfun(@(l)l'<0,lr','UniformOutput',false),mudz,'UniformOutput',false);
end

function [a] = getUpperActives(mudz)
a = cellfun(@(lr)cellfun(@(l)l'>0,lr','UniformOutput',false),mudz,'UniformOutput',false);
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
    me = cellfun(@minus,mU,mL,'UniformOutput',false);
end
function me = plusC(mU,mL)
    f = @plusS;
    me = cellfun(f,mU,mL,'UniformOutput',false);
end
function me = plusS(mU,mL)
    me = cellfun(@plus,mU,mL,'UniformOutput',false);
end

function zDims = getZDims(z)
zDims = cellfun(@(zr)cellfun(@numel,zr),z,'UniformOutput',false);
end

function c = countActiveLines(act)
c = cellfun(@countActiveLinesS,act,'UniformOutput',false);
end
function c = countActiveLinesS(actr)
c = sum(cellfun(@sum,actr));
end

function [Aa,ba,Ra] = applyBuildConstraints(sign,Az,bz,act,az)
[ Aa,ba,Ra] = cellfun(@(A,ld,newConsrk,a)buildActiveConstraints(A,ld,newConsrk,sign,'R',a),Az,bz,act,az,'UniformOutput',false);
end

function zact = extractActiveVector(z,act,sign)
    zact = cellfun(@(xr,ar)cellfun(@(xrk,ark)sign*xrk(ark),xr ,ar,'UniformOutput',false),z,act,'UniformOutput',false);
end

function dz = composeStep(az,dzN,xibar)
	dz = cellfun(@(ar,dr)cellfun(@(z,dz)xibar*z+dz,ar,dr,'UniformOutput',false),az,dzN,'UniformOutput',false);
end
function dzN = nullSpaceStep(Az,duC,ss)
	dzN = cellfun(@(Ar,ssr)cellmtimes(Ar,duC,'lowerTriangular',true,'ci',ssr.ci),Az,ss,'UniformOutput',false);
end

function [feasiblez,lowActivez,upActivez,violationz ] = applyCheckConstraintFeasibility(dz,ldz,udz,feasTol,nCons) 
    [feasiblez,lowActivez,upActivez,violationz ] = cellfun(@(d,l,u)checkConstraintFeasibility(d,l,u,'primalFeasTol',feasTol,'first',nCons),dz,ldz,udz,'UniformOutput',false) ;
    violationz = max([cell2mat(violationz);-inf]);
end

function [mubz,to] = extractDuals(dual,to,ncbz,newConsbz,zDims,mubz,first,sign,jobSchedule)
    imMaster = jobSchedule.imMaster;
    ncbz = bringVariablesMPI(ncbz,jobSchedule);
	if imMaster
    ncbz = cell2mat(ncbz);
    toj = to + cumsum(ncbz);
    fromj = [to+1;toj(1:end-1)+1];
    
    dualR = arrayfun(@(fromr,tor)-sign*dual(fromr:tor),fromj,toj,'UniformOutput',false);
    else
	dualR = [];
	end
    dualR = distributeVariablesMPI(dualR,jobSchedule);
    
    %spmd
    d = applyExtractCompressIneq(dualR,newConsbz,zDims);
    %end
    
    if first
        mubz = d;
    else
        %spmd
        mubz = plusC(mubz,d);
        %end
    end
	if imMaster
	to = toj(end);
    else
	to = 0;
	end
end

function d = applyExtractCompressIneq(dualR,newConsbz,zDims)

    d = cellfun(@(dualRr,newConsbzr,zDimsr)extractCompressIneq(dualRr,newConsbzr,zDimsr,true),dualR,newConsbz,zDims,'UniformOutput',false);

end

function [mubz,to] = extractDualsS(dual,to,ncbz,newConsbz,zDims,mubz,first,sign)

    from = to + 1;
    to = to + ncbz;
    d = extractCompressIneq(-sign*dual(from:to),newConsbz,zDims,true);
    if first
        mubz = d;
    else
        mubz = plusS(mubz,d);
    end

end

function lag = gradientLagrangian(mudz,Az,ss)
if isempty(mudz)
    lag = 0;
else
    lag = catAndSum(cellfun(@(mud,A,ssr)cell2mat(cellmtimesT( mud,A ,'lowerTriangular',true,'ci',ssr.ci,'columnVector',false)),mudz,Az,ss,'UniformOutput',false));
end    
end

function err = checkSolutionFeasibility(P)

dxP = P.Model.A*P.Solution.x;
err = max([P.Model.lhs - dxP;dxP - P.Model.rhs;0]);
end