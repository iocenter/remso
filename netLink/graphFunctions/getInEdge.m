function rowEdge = getInEdge(A,j)
%% For a given vertex, this function returns the row identifying the pipeline (edge) sending production to it.
    ik = 1;
    while (A(ik,j) == 0)|| (ik == j)
           ik = ik+1;         
    end
    rowEdge = ik;

end
