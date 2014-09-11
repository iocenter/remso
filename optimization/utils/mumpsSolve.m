function [ SOL ] = mumpsSolve( mat,RHS )



% initialization of a matlab MUMPS structure
id = initmumps;
id.SYM = 0;

% here JOB = -1, the call to MUMPS will initialize C 
% and fortran MUMPS structure
id = dmumps(id);
id.JOB = 6;


id.ICNTL(3) = 0;
id.ICNTL(4) = 1;


%id.ICNTL(7) = 2; % ordering options!

id.RHS = RHS;

%call to mumps
id = dmumps(id,mat);


SOL = id.SOL;

% kill mumps
id.JOB = -2;
id = dmumps(id);


end

