function [  ] = printLogLine(it,nameList,varList,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

opt = struct('logName','convergenceDetail.txt');
opt = merge_options(opt, varargin{:});

if it ==1
    fid = fopen(opt.logName,'w');
else
	fid = fopen(opt.logName,'a');
end

ePlusSpace = 11;

spaceforIt = 3;
digitsK = numel(num2str(it));
spaceIt = repmat(' ',1,spaceforIt-digitsK);


if mod(it,10) == 1
    header = true;
else
    header = false;
end

if header
    fout = ['it '];
    for k = 1:numel(nameList)-1
       fout = [fout nameList{k} repmat(' ',1,ePlusSpace- numel(nameList{k}))]; 
    end
    fout = [fout nameList{end} '\n']; 
    displayMessage( fout,{},'fid',fid,'std',false);
end



fout = ['%i',spaceIt,repmat('%+.3e ',1,numel(varList)-1),'%+.3e\n'];


displayMessage( fout,[it varList],'fid',fid,'std',false);


fclose(fid);
end

