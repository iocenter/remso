function A = runNetwork(A)
%% runs a full simulation for the whole production network
    A = initNetwork(A);
    
    for i = 1:length(A)  % only the diagonal contains vertices
        if A(i,i) == -2 %% producer
            A = runProducer(A,i);                
        elseif A(i,i) == 2 %% injector
            A = runInjector(A,i);
        end        
    end   
end

function A = runProducer(A,i)
%% This function performs a flow simulation for a producer
    ik=i;
    vert = A(ik,ik);
    while (vert.type ~= -1) %% production ending vertex
        
        %% updating flows and pressure loss in the pipeline
        colEdge = getOutEdge(A, ik);
        dp = dpin_uphill(A(ik, colEdge), vert);
        A(ik,colEdge).qoE = vert.qoV;
        A(ik,colEdge).qwE = vert.qwV;
        A(ik,colEdge).qgE = vert.qgV;
        A(ik,colEdge).dpE = dp;        
    
        %% updating flows and pressure of the ending vertex
        A(colEdge,colEdge).qoV =  A(ik,colEdge).qoE;
        A(colEdge,colEdge).qwV =  A(ik,colEdge).qwE;
        A(colEdge,colEdge).qgV =  A(ik,colEdge).qgE;
        A(colEdge,colEdge).pressure = vert.pressure - dp;
        
        ik = colEdge;
        vert = A(ik,ik);
    end 
end

function A = runInjector(A,j)
%% This function perfoms a simulation for an injection well
    jk = j;
    vert = A(jk,jk);
    while (vert.type ~= 1) %% injection ending vertex
        rowEdge = getInEdge(A, jk);
        dp = dpout_downhill(A(rowEdge, jk), vert);
                
        A(rowEdge,jk).qoE = vert.qoV;
        A(rowEdge,jk).qwE = vert.qwV;
        A(rowEdge,jk).qgE = vert.qgV;
        A(rowEdgeik,jk).dpE = dp;        
        
        
        A(rowEdgeik,rowEdgeik).qoV = vert.qoV;
        A(rowEdgeik,rowEdgeik).qwV= vert.qwV;
        A(rowEdgeik,rowEdgeik).qgV = vert.qgV;
        A(rowEdgeik,rowEdgeik).pressure =  vert.pressure - dp;
        
        jk = rowEdge;
        vert = A(jk,jk);
    end
end