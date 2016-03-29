function pipe = newPipeline(varargin)
% Creates a struct to represent a pipeline        

   if nargin < 4        
        error('Missing Parameters to invoke function newPipeline ');
    end

   opt = struct('diam',[],'len',[],'ang',[],'temp',[], 'roughness', 2.8*10^-5*meter);
   pipe = merge_options(opt, varargin{:});
end

