function [qw,qo,bhp,t0,tf] = perfVariables(W, fluid, simRes)



%-----------------------------------------------


if simRes(1).timeInterval*[-1;1] == 0
    simRes = simRes(2:end);
end

numSteps = numel(simRes);

qw = cell(numSteps,1);
qo = cell(numSteps,1);
bhp = cell(numSteps,1);
t0 = zeros(numSteps,1);
tf = zeros(numSteps,1);

mob = cell([1, 1 ]);


wellSys   = [ W.S ];
Dw        = blkdiag( wellSys.D );
perfSign = Dw*vertcat(W.sign);
wellCells = vertcat( W.cells );

injCells = perfSign == 1;
prodCells = perfSign == -1;




for step = 1 : numSteps,
    t0(step) = simRes(step).timeInterval(1);
    tf(step) = simRes(step).timeInterval(2);
    
    resSol  = simRes(step).resSol;
    wellSol = simRes(step).wellSol;

	wellRates = vertcat(wellSol.flux);
    wellSats  = resSol.s( wellCells );

    [mob{:}] = mobilities(struct('s', wellSats), fluid);

    Lt  = sum(mob{1}, 2);
    f   = bsxfun(@rdivide, mob{1}, Lt);

    f_w = f(:,1);
    f_o = 1 - f_w;


    qw{step} = [ wellRates(injCells)    ; -wellRates(prodCells).*f_w(prodCells) ] ;
	qo{step} = [ zeros(sum(injCells),1) ; -wellRates(prodCells).*f_o(prodCells) ] ;

    bhp{step} = vertcat(wellSol.pressure);

end


end

%--------------------------------------------------------------------------

function [wellRates, rateSigns] = getRates(W, wellSol)
   wellRates = vertcat(wellSol.flux);
   wellSys   = [ W.S ];
   Dw        = blkdiag( wellSys.D );
   wellSigns = ones( numel(W), 1 );
   totRates  = Dw'*wellRates;
   wellSigns( totRates < 0 ) = -1;

   for k = 1:numel(W),
      if ~isempty(W(k).sign) % override if sign is given expl
         wellSigns(k) = W(k).sign;
      end
   end

   rateSigns = Dw*wellSigns;
end

%--------------------------------------------------------------------------

function [mob, dmob, dmob2] = mobilities(state, fluid)
   %output/derivatives should be wrt s_w
   mu = fluid.properties(state);
   s  = fluid.saturation(state);
   [kr{1:nargout}] = fluid.relperm(s, state);

   %        \lambda_i in varargout{1}.
   % (d/ds) \lambda_i in varargout{2}.  Returned iff requested.
   %
   mob = bsxfun(@rdivide, kr{1}, mu);
   if nargout > 1
       dmob = bsxfun(@rdivide, kr{2}(:, [1 end]), mu);
       dmob(:, 2) = -dmob(:,2);
   end
    if nargout > 2
       dmob2 = bsxfun(@rdivide, kr{3}(:, [1 end]), mu);
       %dmob2(:, 2) = -dmob(:,2);
   end
   %kr = cellfun(@(x)x(:,[1 end]), kr, 'UniformOutput', false);
   %varargout = cellfun(@(n) bsxfun(@rdivide, n, mu), kr, ...
   %                    'UniformOutput', false);
end
