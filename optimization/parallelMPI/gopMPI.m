function [ mat ] = gopMPI( op,mat,jobSchedule )
% valid op
% op='+'
% op='*' 
% op='M' --> max
% op='N' --> min
% op='&'
% op='|'

matSize = size(mat);
n = numel(mat);
matV = full(reshape(mat,n,1));

matV=NMPI_Reduce(matV,n,op,jobSchedule.Master_rank);

mat = reshape(matV,matSize(1),matSize(2));

end

