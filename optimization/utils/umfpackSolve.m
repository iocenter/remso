function [ x ] = umfpackSolve(J,b)

[L,U,P,Q] = lu(J);
x = Q * (U \ (L \ (P * b))) ;

end

