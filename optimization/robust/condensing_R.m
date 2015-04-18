function varargout = condensing_R(x,u,v,s,ss,varargin)



opt = struct('simVars',[],'uRightSeeds',[],'computeCorrection',false,'computeNullSpace',true,'xd',[],'vd',[],'sd',[],'eta',0.9);
opt = merge_options(opt, varargin{:});

nR = numel(ss);

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
withAlgs = true;


xs = cell(nR,1);
vs = cell(nR,1);


ax = cell(nR,1);
Ax = cell(nR,1);
av = cell(nR,1);
Av = cell(nR,1);



for kr = 1:nR
    [xs{kr},vs{kr},xd{kr},vd{kr},ax{kr},Ax{kr},av{kr},Av{kr}] = condensing(x{kr},u,v{kr},ss{kr},...
        'simVars',simVars{kr},...
        'uRightSeeds',uRightSeeds,...
        'computeCorrection',computeCorrection,...
        'computeNullSpace',computeNullSpace,...
        'xd',xd{kr},...
        'vd',vd{kr},...
        'withAlgs',withAlgs) ;
end



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
        vDims = cellfun(@(z)cellfun(@numel,z),v,'UniformOutput',false);
        xDims = cellfun(@(z)cellfun(@numel,z),x,'UniformOutput',false);
        uDims = cellfun(@numel,u);
        
        seedSizes = [size(opt.uRightSeeds{1},2),1];
        sumSeedSizes = sum(seedSizes);
        vRightSeed = cellfun(@(A,a,Dimsr)mat2cell([cell2mat(A),cell2mat(a)],Dimsr,sumSeedSizes),Av,av,vDims,'UniformOutput',false);
        xRightSeed = cellfun(@(A,a,Dimsr)mat2cell([cell2mat(A),cell2mat(a)],Dimsr,sumSeedSizes),Ax,ax,xDims,'UniformOutput',false);
        uRightSeed = mat2cell([cell2mat(opt.uRightSeeds),zeros(sum(uDims),1)],uDims,sumSeedSizes);
        
        [s2,sJac] = realization2s(x,u,v,ss,'partials',true,'vRightSeed',vRightSeed,'xRightSeed',xRightSeed,'uRightSeed',uRightSeed,'eta',opt.eta);
        
        % recover data
        As = sJac.J(:,1:end-1);
        dsda = sJac.J(:,end);
        
        if ~givenRangeRHS
            sd = s2-s;
        end
        as = sd+dsda;
        
    else
        
        [s2,sJac] = realization2s(x,u,v,ss,'partials',true,'eta',opt.eta);
        
        % correction on s given the correcion on v (av)
        dsdav = cellfun(@(sJvi,avi)cell2mat(sJvi)*cell2mat(avi),sJac.Jv',av,'UniformOutput',false);
        dsdav = sum(cat(3,dsdav{:}),3);
        dsdax = cellfun(@(sJvi,avi)cell2mat(sJvi)*cell2mat(avi),sJac.Jx',ax,'UniformOutput',false);
        dsdax = sum(cat(3,dsdax{:}),3);
        
        % correction on the stocastic values
        if ~givenRangeRHS
            sd = s2-s;
        end
        as = sd+dsdav+dsdax;
        
        % total derivative of the stochastic variables w.r.t the controls
        Asv = cellfun(@(J,A,ssr)...
            cell2mat(cellmtimesT( J,A,'lowerTriangular',true,'ci',ssr.ci,'columnVector',false)),...
            sJac.Jv',Av,ss,...
            'UniformOutput',false);
        Asx = cellfun(@(J,A,ssr)...
            cell2mat(cellmtimesT( J,A,'lowerTriangular',true,'ci',ssr.ci,'columnVector',false)),...
            sJac.Jx',Ax,ss,...
            'UniformOutput',false);
        
        uDims = cellfun(@numel,u);
       
        As = mat2cell(sum(cat(3,Asv{:}),3)+sum(cat(3,Asx{:}),3)+cell2mat(sJac.Ju),numel(s),uDims);
        
    end
    
elseif computeNullSpace && ~computeCorrection
    
    if size(uRightSeeds{1},1)>0
        
        % TODO: verify that eta is always passed on!
        [s2,sJac] = realization2s(x,u,v,ss,'partials',true,'vRightSeed',Av,'xRightSeed',Ax,'uRightSeed',opt.uRightSeeds,'eta',opt.eta);
        
        % recover data
        As = sJac.J;
        
        if ~givenRangeRHS
            sd = s2-s;
        end
        
    else
        
        [s2,sJac] = realization2s(x,v,u,ss,'partials',true,'eta',opt.eta);
        
        if ~givenRangeRHS
            sd = s2-s;
        end
        
        % total derivative of the stochastic variables w.r.t the controls
        Asv = cellfun(@(Jr,Avr,ssr)...
            cell2mat(cellmtimesT( Jr,Avr,'lowerTriangular',true,'ci',ssr.ci,'columnVector',false)),...
            sJac.Jv',Av,ss,...
            'UniformOutput',false);
        Asx = cellfun(@(J,A,ssr)...
            cell2mat(cellmtimesT( J,A,'lowerTriangular',true,'ci',ssr.ci,'columnVector',false)),...
            sJac.Jx',Ax,ss,...
            'UniformOutput',false);
        
        uDims = cellfun(@numel,u);
        
        As = mat2cell(sum(cat(3,Asv{:}),3)+sum(cat(3,Asx{:}),3)+cell2mat(sJac.Ju),numel(s),uDims);
        
    end
    
elseif ~computeNullSpace && computeCorrection 
             uDims = cellfun(@numel,u);
  
    uRightSeed = cellfun(@(ui)zeros(numel(ui),1),u,'UniformOutput',false);
    [s2,sJac] = realization2s(x,u,v,ss,'partials',true,'vRightSeed',av,'xRightSeed',ax,'uRightSeed',uRightSeed,'eta',opt.eta);
    
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