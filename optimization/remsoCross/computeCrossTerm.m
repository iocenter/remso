function [ w,stepY ] = computeCrossTerm(x,u,v,ax,av,gbarZ,ss,obj,mudx,mudu,mudv,lbx,lbv,ubx,ubv,withAlgs,varargin)

opt = struct('xs',[],'vs',[]);
opt = merge_options(opt, varargin{:});

    [stepY] = maximumStepLength([x;v],[ax;av],[lbx;lbv],[ubx;ubv],'tol',0);
    if stepY == 0
        axN = cellfun(@(x)-1*x,ax,'UniformOutput',false);
        if withAlgs
            avN = cellfun(@(x)-1*x,av,'UniformOutput',false);
        end
        [stepY] = maximumStepLength([x;v],[axN;avN],[lbx;lbv],[ubx;ubv],'tol',0);
        stepY = -stepY;
    end
    if stepY == 0
        % Impossible to compute the correction step by finite differences
        w = cellfun(@(xx)zeros([size(xx,2),size(xx,1)]),u','UniformOutput',false);
    else

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
        
        [xsR,vsR,~,convergedR,simVarsR,uslicedR] = simulateSystem(xR,u,ss,'gradients',false,'guessX',xRG,'guessV',vRG);

        if isa(xsR,'Composite')
            xsR = bringVariables(xsR,ss.jobSchedule);
            vsR = bringVariables(vsR,ss.jobSchedule);
        end
        
        xdR = cellfun(@(xsi,xi)xsi-xi,xsR,xR,'UniformOutput',false);
        vdR = cellfun(@(vsi,vi)vsi-vi,vsR,vR,'UniformOutput',false);
        
        
        [fR,objPartialsR] = obj(xR,u,vR,'gradients',true,'usliced',uslicedR);
        
        gbarRdx.Jx =  cellfun(@(Jz,m)(Jz+m'),objPartialsR.Jx,mudx','UniformOutput',false);
        gbarRdx.Ju =  cellfun(@(Jz,m)(Jz+m'),objPartialsR.Ju,mudu','UniformOutput',false);
        if withAlgs
            gbarRdx.Jv = cellfun(@(Jz,m)(Jz+m'),objPartialsR.Jv,mudv','UniformOutput',false);
        end
        
        [~,gradUY,convergedR,~,~,~,~ ] = simulateSystemZ(u,xdR,vdR,ss,[],'gradients',true,'guessV',xR,'guessX',vR,'simVars',simVarsR,'JacTar',gbarRdx);
        
        if stepY == 1
            w = cellfun(@minus,gradUY,gbarZ,'UniformOutput',false);
        else
            w = cellfun(@(x1,x2)(1/stepY)*(x1-x2),gradUY,gbarZ,'UniformOutput',false);
        end
    end
    


end

