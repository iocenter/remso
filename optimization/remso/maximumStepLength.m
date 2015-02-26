function [maxStep,dx] = maximumStepLength(x,dx,lbxH,ubxH,varargin)

opt = struct('tol',1e-5);
opt = merge_options(opt, varargin{:});

%  given   lbxH <= x <= ubxH
%  solve   max  l  
%           l     
%               s.t  lbxH <= x + l*dx <= ubxH 
%                     l <= 1


xDims = cellfun(@(xi)size(xi,1),x);

x = cell2mat(x);
dx = cell2mat(dx);
lbxH = cell2mat(lbxH);
ubxH = cell2mat(ubxH);

% Determine constraints being violated
ubVVal = ubxH-x-dx;
ubV = ubVVal < 0;


% if violating dx is smaller than the tolerance, chopp it to satisfy boundary
dxTolBool = ubVVal(ubV)>-opt.tol;
ubVT = false(size(ubV));
ubVT(ubV) = dxTolBool;
dx(ubVT) = ubxH(ubVT)-x(ubVT);                
ubV(ubVT) = false;

lMaxub = (ubxH(ubV)-x(ubV))./dx(ubV);

lbVVal = lbxH-x-dx;
lbV =  lbVVal > 0 ;

% if violating dx is smaller than the tolerance, choop it to satisfy boundary
dxTolBool =  lbVVal(lbV)<opt.tol;
lbVT = false(size(lbV));
lbVT(lbV) = dxTolBool;
dx(lbVT) = lbxH(lbVT)-x(lbVT);                
lbV(lbVT) = false;

      
lMaxlb = (lbxH(lbV)-x(lbV))./dx(lbV);


lMaxubT =  min([inf;lMaxub]);    
lMaxlbT =  min([inf;lMaxlb]);

maxStep = min([lMaxlbT,lMaxubT,1]);

assert(maxStep >= 0, 'Check incumbent x and hard constraint bounds')

if maxStep == 0
    warning('maxStep == 0, you may want to increase a bit opt.tol. The QP solver is not accurate enough');
end

dx = mat2cell(dx,xDims,1);
   


end

