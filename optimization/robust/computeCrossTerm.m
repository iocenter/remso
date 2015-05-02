function [ w,stepY ] = computeCrossTerm(x,u,v,s,ax,av,as,gbarZ,sss,obj,mudx,mudu,mudv,muds,lbx,lbv,lbs,ubx,ubv,ubs,varargin)

opt = struct('xs',[],'vs',[],'s2',[],'minStep',1e-6,'tol',sqrt(eps));
opt = merge_options(opt, varargin{:});

tol = opt.tol;

stepY = maximumStepLength({s},{as},{lbs},{ubs},'tol',tol,'debug',false);
if stepY > 0
    
        stepYR = maxStepCalc(x,v,ax,av,lbx,lbv,ubx,ubv,tol);

    
    stepY = min(stepYR,stepY);
end


if stepY == 0  % try to compute in the negative side

    asN = -as;
    stepY = maximumStepLength({s},{asN},{lbs},{ubs},'tol',tol,'debug',false);
    if stepY > 0

    	axN = uMinus(ax);
        avN = uMinus(av);
        stepYR = maxStepCalc(x,v,axN,avN,lbx,lbv,ubx,ubv,tol);

        
        stepY = min(stepYR,stepY);
    end
    stepY = -stepY;
end


if stepY == 0
    % Impossible to compute the correction step by finite differences
    w = cellfun(@(xx)zeros([size(xx,2),size(xx,1)]),u','UniformOutput',false);
else
    convergedR = false;
    
    while ~all(convergedR) && abs(stepY) > opt.minStep
        
        
        if stepY == 1
            xR = plusR(x,ax);
            vR = plusR(v,av);
            sR = s+as;
            
            xRG = xR;
            vRG = vR;
            
        else
            xR = plusRI(x,ax,stepY);
            vR = plusRI(v,av,stepY);
            sR = s+as*stepY;

            
            xRG = [];
            vRG = [];
            
            
            if ~isempty(opt.xs) && ~isempty(opt.vs)
                xs = opt.xs;
                vs = opt.vs;               
                
                xRG = buildGuessC(xs,x,ax,stepY);
                
                vRG = buildGuessC(vs,v,av,stepY);
                
                
                
            end
            
        end
        try                  
            [xsR,vsR,s2R,~,convergedR,simVarsR,uslicedR] = simulateSystem_R(xR,u,vR,sss,'gradients',false,'guessX',xRG,'guessV',vRG);
            if ~all(convergedR)
                stepY = stepY/2;
            end
        catch
            warning(['Simulation fails when computing cross-term. Reducing step= ' num2str(stepY)])
            convergedR = false;
            stepY = stepY/2;
        end
        
        
    end
    if all(convergedR)
        
        
        [fR,objPartialsR] = obj(sR,u,'gradients',true);
        
        gbarRdx.Jx =  mudx;
        gbarRdx.Ju =  plusC(objPartialsR.Ju,mudu);
        gbarRdx.Jv = mudv;
        gbarRdx.Js = objPartialsR.Js+cell2mat(muds);

        
        [gradUY] = simulateSystemZ_R(u,xR,vR,sss,gbarRdx,simVarsR);
        
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

function maxStepxv = maxStepCalc(x,v,ax,av,lbx,lbv,ubx,ubv,tol)
    maxStepxv = min([cellfun(@(xr,vr,axr,avr,lbvr,ubvr)...
        maximumStepLength([xr;vr],[axr;avr],[lbx;lbvr],[ubx;ubvr],'tol',tol,'debug',false),...
        x,v,ax,av,lbv,ubv);inf]);
end

function zm = uMinus(z)
    zm = cellfun(@(zr)cellfun(@uminus,zr,'UniformOutput',false),z,'UniformOutput',false);
end


function zm = minusR(z1,z2)
    zm = cellfun(@minusC,z1,z2,'UniformOutput',false);
end
function zm = minusC(z1,z2)
    zm = cellfun(@minus,z1,z2,'UniformOutput',false);
end

function zm = plusR(z1,z2)
    zm = cellfun(@plusC,z1,z2,'UniformOutput',false);
end
function zm = plusC(z1,z2)
    zm = cellfun(@plus,z1,z2,'UniformOutput',false);
end

function zm = plusRI(z1,z2,l)
    f = @(z1r,z2r)plusCI(z1r,z2r,l);
    zm = cellfun(f,z1,z2,'UniformOutput',false);
end
function zm = plusCI(z1,z2,l)
    zm = cellfun(@(z1r,z2r)z1r+z2r*l,z1,z2,'UniformOutput',false);
end


function zRG = buildGuessC(zs,z,az,stepY)
    f= @(zsr,zr,azr)buildGuessM(zsr,zr,azr,stepY);
   zRG = cellfun(f,zs,z,az,'UniformOutput',false);
end
function zRG = buildGuessM(zs,z,az,stepY)
   zRG = cellfun(@(x1,x2,x3)x1*(1-stepY)+(x2+x3)*stepY,zs,z,az,'UniformOutput',false);
end
