function [ x ] = umfpackSolve(J,b)

[L,U,P,Q,R] = lu(J);
x = Q * (U \ (L \ (P * (R\b)))) ;

end

