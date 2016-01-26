function colEdge = getOutEdge(A,i)
%% For a given vertex, this function returns the column identifying the pipeline (edge) flowing its production.
    jk = 1;
    while (A(i,jk) ==0) || (i==jk)
        jk = jk+1;
    end
    colEdge = jk;
end