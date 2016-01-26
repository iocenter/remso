function pipe = newPipeline(varargin)
% Creates a struct to represent a pipeline        

   if nargin < 4        
        error('Missing Parameters to invoke function newPipeline ');
    end

   opt = struct('diam',[],'len',[],'ang',[],'temp',[]);
   pipe = merge_options(opt, varargin{:});
end

