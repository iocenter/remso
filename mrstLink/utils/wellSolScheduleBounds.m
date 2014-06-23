function [wellSolMax,wellSolMin] = wellSolScheduleBounds(wellSol,varargin)
%
% fill a wellSol mock object with maximums and minimums according to the
% input parameters
%

opt     = struct('maxProd',struct('WRAT',inf,'ORAT',inf,'GRAT',inf,'BHP',inf),...
                 'minProd',struct('WRAT',0,  'ORAT',0,  'GRAT',0,  'BHP',0),...
                 'maxInj', struct('WRAT',inf,'ORAT',inf,'GRAT',inf,'BHP',inf),...
                 'minInj', struct('WRAT',0,  'ORAT',0,  'GRAT',0,  'BHP',0),...
    'useWellLims',[]);

opt = merge_options(opt, varargin{:});


if ~isempty(opt.useWellLims)
	%TODO: 
	error('Not Implemented');
end


wellSolMax = wellSol;
wellSolMin = wellSol;

for w = 1:numel(wellSol)
    switch wellSol(w).sign
        case -1
            wellSolMax(w).bhp =  opt.maxProd.BHP;
            wellSolMax(w).qOs = -opt.minProd.ORAT;
            wellSolMax(w).qWs = -opt.minProd.WRAT;
            wellSolMax(w).qGs = -opt.minProd.GRAT;
            
            wellSolMin(w).bhp =  opt.minProd.BHP;
            wellSolMin(w).qOs = -opt.maxProd.ORAT;
            wellSolMin(w).qWs = -opt.maxProd.WRAT;
            wellSolMin(w).qGs = -opt.maxProd.GRAT;
        case 1
            wellSolMax(w).bhp = opt.maxInj.BHP;
            wellSolMax(w).qOs = opt.maxInj.ORAT;
            wellSolMax(w).qWs = opt.maxInj.WRAT;
            wellSolMax(w).qGs = opt.maxInj.GRAT;
           
            
            wellSolMin(w).bhp = opt.minInj.BHP;
            wellSolMin(w).qOs = opt.minInj.ORAT;
            wellSolMin(w).qWs = opt.minInj.WRAT;
            wellSolMin(w).qGs = opt.minInj.GRAT;

        otherwise
            error('what');
    end
    
end




end