function [ hhxvu ] = buildFullHessian( x,v,u,obj,ss,lambdaX,lambdaV,varargin )


opt = struct('pert',1e-5);
opt = merge_options(opt, varargin{:});

% the values of mu and bounds do not affect the Hessian.

muU = cellfun(@(zi)rand(size(zi))',u','UniformOutput',false);
muX = cellfun(@(zi)rand(size(zi))',x','UniformOutput',false);
muV = cellfun(@(zi)rand(size(zi))',v','UniformOutput',false);


lbx = cellfun(@(zi)min(zi+rand(size(zi))-0.5,zi),x,'UniformOutput',false);
ubx = cellfun(@(zi,lbxi)max(max(zi+rand(size(zi))-0.5,lbxi+0.1),zi),x,lbx,'UniformOutput',false);

lbv = cellfun(@(zi)min(zi+rand(size(zi))-0.5,zi),v,'UniformOutput',false);
ubv = cellfun(@(zi,lbxi)max(max(zi+rand(size(zi))-0.5,lbxi+0.1),zi),v,lbv,'UniformOutput',false);

lbu = cellfun(@(zi)min(zi+rand(size(zi))-0.5,zi),u,'UniformOutput',false);
ubu = cellfun(@(zi,lbxi)max(max(zi+rand(size(zi))-0.5,lbxi+0.1),zi),u,lbu,'UniformOutput',false);




[ lagF,lagG] = lagrangianF( u,x,v,lambdaX,lambdaV,muU,muX,muV,obj,ss,lbx,lbv,lbu,ubx,ubv,ubu,'gradients',true);



% construct hessian
pert = opt.pert;

hGz = cell(numel(x),1);

fields = fieldnames(lagG);

parfor zci = 1:numel(x)
    
    
    hGzi = cell(numel(x{zci}),1); 
    for zcii = 1:numel(x{zci})
        
        zi = x;
        zi{zci}(zcii) = zi{zci}(zcii) + pert;
        
        [~,hGzi{zcii}] = lagrangianF( u,zi,v,lambdaX,lambdaV,muU,muX,muV,obj,ss,lbx,lbv,lbu,ubx,ubv,ubu,'gradients',true);
        
    end
    
    jj = [hGzi{:}];
    
    for k = 1:numel(fields)
   
        hGz{zci}.(fields{k}) = cell2mat(cellfun(@(xx)cell2mat(xx),{jj.(fields{k})}','UniformOutput',false));

    end
    
end

hhxvu = cell(3,3);

jj = [hGz{:}];
for k = 1:numel(fields)
    hhz.(fields{k}) = cell2mat({jj.(fields{k})}');
    hhz.(fields{k}) = (hhz.(fields{k}) - ones(size(hhz.(fields{k}),1),1)*cell2mat(lagG.(fields{k})))/pert;
end

hhxvu{1,1} = hhz.Jx;
hhxvu{1,2} = hhz.Jv;
hhxvu{1,3} = hhz.Ju;

hGz = cell(numel(v),1);

parfor zci = 1:numel(v)
    
    
    hGzi = cell(numel(v{zci}),1); 
    for zcii = 1:numel(v{zci})
        
        zi = v;
        zi{zci}(zcii) = zi{zci}(zcii) + pert;
        
        [~,hGzi{zcii}] = lagrangianF( u,x,zi,lambdaX,lambdaV,muU,muX,muV,obj,ss,lbx,lbv,lbu,ubx,ubv,ubu,'gradients',true);
        
    end
    
    jj = [hGzi{:}];
    
    for k = 1:numel(fields)
   
        hGz{zci}.(fields{k}) = cell2mat(cellfun(@(xx)cell2mat(xx),{jj.(fields{k})}','UniformOutput',false));

    end
    
end


jj = [hGz{:}];
for k = 1:numel(fields)
    hhz.(fields{k}) = cell2mat({jj.(fields{k})}');
    hhz.(fields{k}) = (hhz.(fields{k}) - ones(size(hhz.(fields{k}),1),1)*cell2mat(lagG.(fields{k})))/pert;
end

hhxvu{2,1} = hhz.Jx;
hhxvu{2,2} = hhz.Jv;
hhxvu{2,3} = hhz.Ju;


hGz = cell(numel(u),1);

parfor zci = 1:numel(u)
    
    
    hGzi = cell(numel(u{zci}),1); 
    for zcii = 1:numel(u{zci})
        
        zi = u;
        zi{zci}(zcii) = zi{zci}(zcii) + pert;
        
        [~,hGzi{zcii}] = lagrangianF( zi,x,v,lambdaX,lambdaV,muU,muX,muV,obj,ss,lbx,lbv,lbu,ubx,ubv,ubu,'gradients',true);
        
    end
    
    jj = [hGzi{:}];
    
    for k = 1:numel(fields)
   
        hGz{zci}.(fields{k}) = cell2mat(cellfun(@(xx)cell2mat(xx),{jj.(fields{k})}','UniformOutput',false));

    end
    
end


jj = [hGz{:}];
for k = 1:numel(fields)
    hhz.(fields{k}) = cell2mat({jj.(fields{k})}');
    hhz.(fields{k}) = (hhz.(fields{k}) - ones(size(hhz.(fields{k}),1),1)*cell2mat(lagG.(fields{k})))/pert;
end

hhxvu{3,1} = hhz.Jx;
hhxvu{3,2} = hhz.Jv;
hhxvu{3,3} = hhz.Ju;


uDims = cellfun(@(ui)numel(ui),u);
xDims = cellfun(@(ui)numel(ui),x);
vDims = cellfun(@(ui)numel(ui),v);


hhxvu = cell2mat(hhxvu);
hhxvu = (hhxvu + hhxvu')/2;  %% make sure the approximation is symmetric

dims = [sum(xDims) sum(vDims) sum(uDims)];

hhxvu = mat2cell(hhxvu,dims',dims);



%{
% check hessian with adimat

uDims = cellfun(@(ui)numel(ui),u);
xDims = cellfun(@(ui)numel(ui),x);
vDims = cellfun(@(ui)numel(ui),v);

xIndex = 1:sum(xDims);
vIndex = sum(xDims)+1:sum(xDims)+sum(vDims);
uIndex = sum(xDims)+sum(vDims)+1:sum(xDims)+sum(vDims)+sum(uDims);

lagFun =@(z) lagrangianF( mat2cell(z(uIndex),uDims,1),...
    mat2cell(z(xIndex),xDims,1),...
    mat2cell(z(vIndex),vDims,1),...
    lambdaX,lambdaV,muU,muX,muV,obj,ss,lbx,lbv,lbu,ubx,ubv,ubu,'gradients',false);



Hv = admDiffFD(@admDiffFD, 1, lagFun, 1, cell2mat([x;v;u]), admOptions('i', [3], 'fdStep', 10*sqrt(sqrt(eps))))

norm(Hv-cell2mat(hhxvu))

%}


end

