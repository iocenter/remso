function varargout = condensing_R(x,u,v,s,sss,varargin)



opt = struct('simVars',[],'uRightSeeds',[],'computeCorrection',false,'computeNullSpace',true,'xd',[],'vd',[],'sd',[]);
opt = merge_options(opt, varargin{:});

ss = sss.ss;
nR = sss.nR;

simVars = opt.simVars;
if isempty(simVars)
    simVars = cell(nR,1);
end
uRightSeeds = opt.uRightSeeds;
computeCorrection = opt.computeCorrection;
computeNullSpace = opt.computeNullSpace;
xd = opt.xd;
if isempty(xd)
    xd = cell(nR,1);
end
vd = opt.vd;
if isempty(vd)
    vd = cell(nR,1);
end



[xs,vs,xd,vd,ax,Ax,av,Av] = applyCondensing(x,u,v,ss,simVars,uRightSeeds,computeCorrection,computeNullSpace,xd,vd);




if isempty(opt.sd)
    sd = [];
    givenRangeRHS = false;
else
    sd = opt.sd;
    givenRangeRHS = true;
end
as = [];
As = [];
if computeNullSpace && computeCorrection
    if ~isempty(uRightSeeds)
        
        % merge seeds
        vDims = getDims(v);
        xDims = getDims(x);     
        uDims = getDimsS(u);
        
        seedSizes = [size(opt.uRightSeeds{1},2),1];
        sumSeedSizes = sum(seedSizes);
        xRightSeed = generateSeeds(Ax,ax,xDims,sumSeedSizes);
        vRightSeed = generateSeeds(Av,av,vDims,sumSeedSizes);
        uRightSeed = mat2cell([cell2mat(opt.uRightSeeds),zeros(sum(uDims),1)],uDims,sumSeedSizes);
        
        [s2,sJac] = realization2s(x,u,v,sss,'partials',true,'vRightSeed',vRightSeed,'xRightSeed',xRightSeed,'uRightSeed',uRightSeed);
        
        % recover data
        As = sJac.J(:,1:end-1);
        dsda = sJac.J(:,end);
        
        if ~givenRangeRHS
            sd = s2-s;
        end
        as = sd+dsda;
        
    else
        
        [s2,sJac] = realization2s(x,u,v,sss,'partials',true);
        
        % correction on s given the correcion on v (av)
        Jv = sJac.Jv;
        Jx = sJac.Jx;
        

        dsdax = sJzTimesdz(Jx,ax);
        dsdav = sJzTimesdz(Jv,av);
       
        
        
        % correction on the stocastic values
        if ~givenRangeRHS
            sd = s2-s;
        end
        as = sd+dsdav+dsdax;
        
        % total derivative of the stochastic variables w.r.t the controls
        Jv = sJac.Jv;
        Jx = sJac.Jx;
            Asx = yTimesA(Jx,Ax,ss);
            Asv = yTimesA(Jv,Av,ss);
                     
            Asv = catAndSum(Asv);

            
            Asx = catAndSum(Asx);


        
        uDims = cellfun(@numel,u);
       
        As = mat2cell(Asx+Asv+sJac.Ju,numel(s),uDims);
        
    end
    
elseif computeNullSpace && ~computeCorrection
    
    if ~isempty(uRightSeeds)
        
        [s2,sJac] = realization2s(x,u,v,sss,'partials',true,'vRightSeed',Av,'xRightSeed',Ax,'uRightSeed',uRightSeeds);
        
        % recover data
        As = sJac.J;
        
        if ~givenRangeRHS
            sd = s2-s;
        end
        
    else
        
        [s2,sJac] = realization2s(x,u,v,sss,'partials',true);
        
        if ~givenRangeRHS
            sd = s2-s;
        end
        
       
        % total derivative of the stochastic variables w.r.t the controls
        Jv = sJac.Jv;
        Jx = sJac.Jx;

            Asv = yTimesA(Jv,Av,ss);
         
            Asx = yTimesA(Jx,Ax,ss);
            
            Asv = catAndSum(Asv);

            
            Asx = catAndSum(Asx);



        
        uDims = cellfun(@numel,u);
        
        As = mat2cell(Asv+Asx+sJac.Ju,numel(s),uDims);
        
    end
    
elseif ~computeNullSpace && computeCorrection 
	uDims = cellfun(@numel,u);
  
    uRightSeed = cellfun(@(ui)zeros(numel(ui),1),u,'UniformOutput',false);
    [s2,sJac] = realization2s(x,u,v,sss,'partials',true,'vRightSeed',av,'xRightSeed',ax,'uRightSeed',uRightSeed);
    
    % recover data
    dsdav = sJac.J;    
    if ~givenRangeRHS
        sd = s2-s;
    end
    as = sd+dsdav;
         
else
    error('Unnecessary call to condensing_R')
end




varargout{1} = xs;
varargout{2} = vs;
varargout{3} = s2;
varargout{4} = xd;
varargout{5} = vd;
varargout{6} = sd;
varargout{7} = ax;
varargout{8} = Ax;
varargout{9} = av;
varargout{10} = Av;
varargout{11} = as;
varargout{12} = As;

end
function [xs,vs,xd,vd,ax,Ax,av,Av] = applyCondensing(x,u,v,ss,simVars,uRightSeeds,computeCorrection,computeNullSpace,xd,vd)


nr = numel(ss);

printCounter= true;
fid = 1;


xs = cell(nr,1);
vs = cell(nr,1);
ax = cell(nr,1);
Ax = cell(nr,1);
av = cell(nr,1);
Av = cell(nr,1);

for r = 1:nr
    if printCounter
        printRef = sprintf('%d/%d',r,nr);
    end
    [xs{r},vs{r},xd{r},vd{r},ax{r},Ax{r},av{r},Av{r}] = ...
        condensing(x{r},u,v{r},ss{r},...
        'simVars',simVars{r},...
        'uRightSeeds',uRightSeeds,...
        'computeCorrection',computeCorrection,...
        'computeNullSpace',computeNullSpace,...
        'xd',xd{r},...
        'vd',vd{r},...
        'withAlgs',true,...
        'printCounter',printCounter,...
        'fid',fid,...
        'printRef',printRef);
end

end

function [zDims] = getDims(z)
    zDims = cellfun(@getDimsS,z,'UniformOutput',false);
end
function [zDims] = getDimsS(z)
    zDims = cellfun(@numel,z);
end


function yA = yTimesA(y,A,ss)

yA = cellfun(@(...
    yr,Ar,ssr)...
    cell2mat(cellmtimesT( yr,Ar,'lowerTriangular',true,'ci',ssr.ci,'columnVector',false)),...
    y , A,ss,'UniformOutput',false);

end

function out = catAndSum(M)
if ~isempty(M)
    if iscell(M{1})
        M = cellfun(@cell2mat,M,'UniformOutput',false);
    end
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

else
    out = 0;
end
end

function sJzTdz = sJzTimesdz(sJz,dz)

sJzTdz = cellfun(@(sJvi,avi)cell2mat(sJvi)*cell2mat(avi),sJz,dz,'UniformOutput',false);
sJzTdz = catAndSum(sJzTdz);

end


function seedsZ = generateSeeds(Az,az,zDims,sumSeedSizes)

seedsZ = cellfun(@(A,a,Dimsr)mat2cell([cell2mat(A),cell2mat(a)],Dimsr,sumSeedSizes),Az,az,zDims,'UniformOutput',false);

end
