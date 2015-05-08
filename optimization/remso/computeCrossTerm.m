function [ w,stepY ] = computeCrossTerm(x,u,v,ax,av,gbarZ,ss,obj,mudx,mudu,mudv,lbx,lbv,ubx,ubv,withAlgs,varargin)

opt = struct('xs',[],'vs',[],'minStep',1e-6,'tol',sqrt(eps),'maxIts',2,'backTrack',0.1);
opt = merge_options(opt, varargin{:});

[stepY] = maximumStepLength([x;v],[ax;av],[lbx;lbv],[ubx;ubv],'tol',opt.tol,'debug',false);
if stepY == 0
    axN = cellfun(@(x)-1*x,ax,'UniformOutput',false);
    if withAlgs
        avN = cellfun(@(x)-1*x,av,'UniformOutput',false);
        [stepY] = maximumStepLength([x;v],[axN;avN],[lbx;lbv],[ubx;ubv],'tol',opt.tol,'debug',false);
    else
        [stepY] = maximumStepLength(x,axN,lbx,ubx,'tol',opt.tol,'debug',false);
    end
    stepY = -stepY;
end
if stepY == 0
    % Impossible to compute the correction step by finite differences
    w = cellfun(@(xx)zeros([size(xx,2),size(xx,1)]),u','UniformOutput',false);
else
    convergedR = false;
    its = 0;
    while ~all(convergedR) && abs(stepY) > opt.minStep && its <= opt.maxIts
        if stepY == 1
            xR = cellfun(@(z,dz)z+dz,x,ax,'UniformOutput',false);
            vR = cellfun(@(z,dz)z+dz,v,av,'UniformOutput',false);
            
            xRG = xR;
            vRG = vR;
            
        else
            xR = cellfun(@(z,dz)z+stepY*dz,x,ax,'UniformOutput',false);
            vR = cellfun(@(z,dz)z+stepY*dz,v,av,'UniformOutput',false);
            
            xRG = [];
            vRG = [];
            if ~isempty(opt.xs)
                xRG = cellfun(@(x1,x2,x3)x1*(1-stepY)+(x2+x3)*stepY,opt.xs,x,ax,'UniformOutput',false);
            end
            
            if ~isempty(opt.vs)
                vRG = cellfun(@(x1,x2,x3)x1*(1-stepY)+(x2+x3)*stepY,opt.vs,v,av,'UniformOutput',false);
            end
            
        end
        try
            [xsR,vsR,~,convergedR,simVarsR,uslicedR] = simulateSystem(xR,u,ss,'gradients',false,'guessX',xRG,'guessV',vRG,'withAlgs',withAlgs);
            if ~all(convergedR)
                stepY = stepY*opt.backTrack;
            end
        catch
            warning(['Simulation fails when computing cross-term. Reducing step= ' num2str(stepY)])
            convergedR = false;
            stepY = stepY*opt.backTrack;
        end
        its = its + 1;
    end
    if all(convergedR)
        
        
        [fR,objPartialsR] = obj(xR,u,vR,'gradients',true,'usliced',uslicedR);
        
        gbarRdx.Jx =  cellfun(@(Jz,m)(Jz+m'),objPartialsR.Jx,mudx','UniformOutput',false);
        gbarRdx.Ju =  cellfun(@(Jz,m)(Jz+m'),objPartialsR.Ju,mudu','UniformOutput',false);
        if withAlgs
            gbarRdx.Jv = cellfun(@(Jz,m)(Jz+m'),objPartialsR.Jv,mudv','UniformOutput',false);
        end
        
        [~,gradUY] = simulateSystemZ(u,xR,vR,ss,obj,'simVars',simVarsR,'JacTar',gbarRdx,'withAlgs',withAlgs);

        
        if stepY == 1
            w = cellfun(@minus,gradUY,gbarZ,'UniformOutput',false);
        else
            w = cellfun(@(x1,x2)(1/stepY)*(x1-x2),gradUY,gbarZ,'UniformOutput',false);
        end
    else
        % This is very bad, lets see if we can continue without the
        % cross term
        w = cellfun(@(xx)zeros([size(xx,2),size(xx,1)]),u','UniformOutput',false);
        stepY = 0;
    end
end



end

