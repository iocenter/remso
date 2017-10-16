function [ schedule,Jac ] = controls2Schedule( u,schedule,varargin)
%
%  set the controls in the schedule.  Scale the variables
%
%

opt = struct('uScale',[],'partials',false,'uRightSeeds',[], 'fixedWells', []);
opt = merge_options(opt, varargin{:});

nu = numel(u);

if ~isempty(opt.uScale)
    u = u.*opt.uScale;
end

if ~isempty(opt.fixedWells)      
    nW = numel(schedule.control(1).W);
    nC = numel(schedule.control);
    nFull = nW*nC;
    
    assert(all(arrayfun(@(c) numel(c.W), schedule.control)==nW)); %% all controls have the same number of wells 
    
    fixedWellsCells = cell2mat(arrayfun(@(x) x + opt.fixedWells,  (0:nC-1)'*nW,'UniformOutput',false ));
    controlWells = setdiff(1:nFull, fixedWellsCells);
    
    uFixed = schedule2Controls(schedule);    
    
    uFull = zeros(nFull,1);    
    uFull(fixedWellsCells) = uFixed(fixedWellsCells);    
    uFull(controlWells) = u;
else
   uFull = u;    
   controlWells = 1:nu;
end

[vals,nC,nW] = controls2CellControls(uFull,schedule);

schedule = cellControls2Schedule(vals,schedule);


if opt.partials    
    if ~isempty(opt.uScale)
        Jac = sparse(controlWells,1:nu,opt.uScale, nW*nC, nu);
    else
%         Jac = speye(nu);
         Jac = speye(controlWells,1:nu,1, nW*nC, nu);
    end
    if size(opt.uRightSeeds,1) ~= 0
        Jac = Jac*opt.uRightSeeds;
    end
    Jac = mat2cell(Jac,repmat(nW,nC,1),size(Jac,2));
else
    Jac = [];
end

end


