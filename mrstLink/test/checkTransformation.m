function [e] = checkTransformation(mrstState,remsoState,xScale,activeComponents,fluid,system )
%
%  TODO: This can be for general Transformation functions, and not only for
%  stateMrst2statePsWrGH
%


e = 0;

if ~isempty(mrstState)
    
    
    [ stateVector,invJac] = stateMrst2stateVector(mrstState,...
        'xScale',xScale,...
        'activeComponents',activeComponents,...
        'fluid',fluid,...
        'system',system,...
        'partials',true);
    
    
    [ mrstState1,Jac ] = stateVector2stateMrst( stateVector,...
        'xScale',xScale,...
        'activeComponents',activeComponents,...
        'fluid',fluid,...
        'system',system,...
        'partials',true);
    
    
    ep = norm(mrstState.pressure-mrstState1.pressure);
    es = norm(mrstState.s-mrstState1.s);
    if isfield(mrstState,'rs')
        ers = norm(mrstState.rs-mrstState1.rs);
        erv = norm(mrstState.rv-mrstState1.rv);
    else
        ers = 0;
        erv = 0;
    end
    e = max([e,ep,es,ers,erv]);
    
    
    
    eJ1 = norm(invJac*Jac - eye(size(Jac)));
    eJ2 = norm(Jac*invJac - eye(size(Jac)));
    
    
    e= max([e,eJ1,eJ2]);
    
end
if ~isempty(remsoState)
    
    [ mrstState2,Jac ] = stateVector2stateMrst( remsoState,...
        'xScale',xScale,...
        'activeComponents',activeComponents,...
        'fluid',fluid,...
        'system',system,...
        'partials',true);
    
    [ remsoState2,invJac] = stateMrst2stateVector(mrstState2,...
        'xScale',xScale,...
        'activeComponents',activeComponents,...
        'fluid',fluid,...
        'system',system,...
        'partials',true);
    
    ep2 = norm(remsoState-remsoState2);
    
    e = max([e,ep2]);
    
    
    eJ12 = norm(Jac*invJac - eye(size(Jac)));
    eJ22 = norm(invJac*Jac - eye(size(Jac)));
    
    e= max([e,eJ12,eJ22]);
    
    
end


end

