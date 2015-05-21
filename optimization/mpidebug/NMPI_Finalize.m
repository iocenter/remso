function [ ] = NMPI_Finalize( )

try
   delete('*MPI.ack');
   delete('*MPI.mat');
catch ex
	msgString = getReport(ex);
	display(msgString);  
end

end

