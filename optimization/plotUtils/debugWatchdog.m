function [ ] = debugWatchdog( it,stepType,l,merit,meritGrad,meritDebugInfo,varargin)

opt = struct('logName','logLineSearch.txt','header',true);
opt = merge_options(opt, varargin{:});

if it ==1
    fid = fopen(opt.logName,'w');
else
	fid = fopen(opt.logName,'a');
end

ePlusSpace = 11;

spaceforIt = 4;
digitsK = numel(num2str(it));
spaceIt = repmat(' ',1,spaceforIt-digitsK);

if mod(it,10) == 1
    header = true;
else
    header = false;
end

nameList = {'l','meritVal','meritGrad','armijo','obj','|eq|_1','rho'};

if header && opt.header
    fout = ['it  T '];
    for k = 1:numel(nameList)
       fout = [fout nameList{k} repmat(' ',1,ePlusSpace- numel(nameList{k}))]; 
    end
    fout = [fout '\n']; 
    displayMessage( fout,{},'fid',fid,'std',false);
end


varList = {it,...
          stepType,...
          l(1),...
          merit(1),...
          meritGrad(1),...
          meritDebugInfo{1}.armijoVal,...
          meritDebugInfo{1}.f,...
          meritDebugInfo{1}.eqNorm1,...
          meritDebugInfo{1}.rho};
fout = ['%i',spaceIt,'%s ',repmat('%+.3e ',1,numel(varList)-2),'\n'];
displayMessage( fout,varList,'fid',fid,'std',false);


for k = 2:numel(l)
    varList =     {l(k)-l(1),...
                   merit(k)-merit(1),...
                   meritGrad(k),....
                   meritDebugInfo{k}.armijoVal-meritDebugInfo{1}.armijoVal,...
                   meritDebugInfo{k}.f-meritDebugInfo{1}.f,...
                   meritDebugInfo{k}.eqNorm1-meritDebugInfo{1}.eqNorm1};
           
    fout = [repmat(' ',1,spaceforIt+2),repmat('%+.3e ',1,numel(varList)),'\n'];
    displayMessage( fout,varList,'fid',fid,'std',false);
    
end


fclose(fid);
end

