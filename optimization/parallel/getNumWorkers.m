function nWorkers = getNumWorkers()

if 2 == exist('gcp','file')
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        nWorkers = 0;
    else
        nWorkers = p.NumWorkers;
    end
else
    nWorkers = matlabpool('size');
end



end