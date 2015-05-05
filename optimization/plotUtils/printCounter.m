function [tprinted,iprinted] =  printCounter(from, to, i,label,t0,it0,fid)

if nargin >=5
    dt = toc(t0);
    if dt > 1 && (dt)/(i-it0)*(to-it0)-dt > 1; % avoid bug
        tprinted  = tic; 
        iprinted  = i;
    elseif i == to
        tprinted  = tic; 
        iprinted  = i;        
    else
        tprinted  = t0; 
        iprinted  = it0;
        return
    end

else
    tprinted  = []; 
    iprinted  = [];
    it0 = from;
end
if nargin<7
    fid = 1;
end;




digits = ceil(log10(max(to,from) + 1));
% print first part if i is at the beguining or if it0 is out of interval
% and never when i is at the end
if (i == from ||  ~xor(it0<=from,it0<=to)) && (i ~= to)
    
    tekst = [label ' ' num2str(max(to,from)) ':%0' num2str(digits) '.0f'];
    fprintf(fid,tekst,i);
    % at the end delete everything if someting was printed before (it0 in the
    % range)
elseif i == to 
    if xor(it0<=from,it0<=to)
        if fid == 1
            delete = repmat('\b', 1, numel(label)+digits*2+2);
            fprintf(fid,delete);
        else
            nChar = (numel(label)+digits*2+2);
            fseek(fid, -nChar , 'cof');
            fprintf(fid,repmat(' ',1,nChar));
            fseek(fid, -nChar , 'cof');
        end
    end
    % remove last number and print the current
else
    if fid == 1
        bstext = repmat('\b', 1, digits);
        tekst = [bstext '%0' num2str(digits) '.0f'];
        fprintf(fid,tekst,i);
    else
        fseek(fid, -digits , 'cof');
        tekst = ['%0' num2str(digits) '.0f'];
        fprintf(fid,tekst,i);       
    end
end




