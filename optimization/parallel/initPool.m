function [ ] = initPool( varargin )



if 2 == exist('gcp','file')
    initPool2014a( varargin{:} );
else
    initPool2013b( varargin{:} );
end
    
end

