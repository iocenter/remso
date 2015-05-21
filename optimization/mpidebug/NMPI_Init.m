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

try
   delete('*MPI.ack');
   delete('*MPI.mat');
catch ex
	msgString = getReport(ex);
	display(msgString);  
end

end

