function maxStep = maximumStepLength(x,dx,lbxH,ubxH)
%  given   lbxH <= x <= ubxH
%  solve   max  l  
%           l     
%               s.t  lbxH <= x + l*dx <= ubxH 
%                     l <= 1

ubV = cellfun(@(xi,ubxHi,dxi)ubxHi-xi-dxi < 0 ,x,ubxH,dx,'UniformOutput',false);

lMaxub = cellfun(@(xi,ubxHi,dxi,ubVi)(ubxHi(ubVi)-xi(ubVi))./dxi(ubVi),x,ubxH,dx,ubV,'UniformOutput',false);


assert( all(cellfun(@(xi)all(isfinite(xi)),lMaxub)),...
       'Check  x <= ubH');

lbV = cellfun(@(xi,lbxHi,dxi)lbxHi-xi-dxi > 0 ,x,lbxH,dx,'UniformOutput',false);
   
   
lMaxlb = cellfun(@(xi,lbxHi,dxi,lbVi)(lbxHi(lbVi)-xi(lbVi))./dxi(lbVi),x,lbxH,dx,lbV,'UniformOutput',false);


assert( all(cellfun(@(xi)all(isfinite(xi)),lMaxlb)),...
       'Check lbH <=x');

lMaxubT = min(cellfun(@(xi) min([inf;xi]),lMaxub));    
lMaxlbT = min(cellfun(@(xi) min([inf;xi]),lMaxlb));

maxStep = min([lMaxlbT,lMaxubT,1]);

assert(maxStep > 0, 'Check incumbent x and hard constraint bounds')
   


end

