function [] = NMPI_Init(varargin )
%NMPI_INIT Summary of this function goes here
%   Detailed explanation goes here
global Master_Rank
global num_ranks
global my_rank

Master_Rank = 0;
num_ranks = 1;
my_rank = 0;

if nargin >= 1
    num_ranks = varargin{1};
end
if nargin >= 2
    my_rank = varargin{2};
end

if my_rank == Master_Rank
    restart = false;
    fmat=dir('*_*_MPI.mat');
    if numel(fmat) > 0
        delete('*_*_MPI.mat');
        restart = true;
    end
    fack=dir('*_*_MPI.ack');
    if numel(fack) > 0
        delete('*_*_MPI.ack');
        restart = true;
    end
    if restart
       error('The communication channel was not clean at start.  Now it must be, restart all ranks') 
    end
end

NMPI_Bcast(0,1,Master_Rank,my_rank);


end

