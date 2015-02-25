function [e] = checkTransformation(mrstState,remsoState,toMRST,toRemso )


e = -inf;

if ~isempty(mrstState)
    
    
    [ x1,x2,x3] = toRemso(mrstState);
    
	stateVector = [x1;x2;x3];
               
    invJac = cell2mat(stateVector.jac);
    stateVector = double(stateVector);
    
    [ mrstState1,Jac ] = toMRST(stateVector);
    Jac = cell2mat(Jac);
    
    ep = norm(mrstState.pressure-mrstState1.pressure);
    es = norm(mrstState.s-mrstState1.s);
    if isfield(mrstState,'rs')
        ers = norm(mrstState.rs-mrstState1.rs);
    else
        ers = 0;
    end
    if isfield(mrstState,'rv')
        erv = norm(mrstState.rv-mrstState1.rv);
    else
        erv = 0;
    end
    e = max([e,ep,es,ers,erv]);
    
    
    
    eJ1 = norm(invJac*Jac - eye(size(Jac)));
    eJ2 = norm(Jac*invJac - eye(size(Jac)));
    
    
    e= max([e,eJ1,eJ2]);
    
end
if ~isempty(remsoState)
    
    [ mrstState2,Jac ] = toMRST(remsoState);
	Jac = cell2mat(Jac);

    
    [ x1,x2,x3] = toRemso(mrstState2);
        
	remsoState2 = [x1;x2;x3];
               
	invJac = cell2mat(remsoState2.jac);
    remsoState2 = double(remsoState2);    
       
    
    ep2 = norm(remsoState-remsoState2);
    
    e = max([e,ep2]);
    
    
    eJ12 = norm(Jac*invJac - eye(size(Jac)));
    eJ22 = norm(invJac*Jac - eye(size(Jac)));
    
    e= max([e,eJ12,eJ22]);
    
    
end


end

