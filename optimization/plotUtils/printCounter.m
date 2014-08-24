function [tprinted,iprinted] =  printCounter(from, to, i,label,t0,it0)

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



if to == -1
    digits = 8;
else
    digits = ceil(log10(to + 1));
end
if i == from || (it0 < from && i ~= to)
    if to == -1
        tekst = [label ':%0' num2str(digits) '.0f'];
    else
        tekst = [label  num2str(to) ':%0' num2str(digits) '.0f'];
    end
    fprintf(tekst,i);
elseif i == to
    if it0 >= from
        delete = repmat('\b', 1, numel(label)+digits*2+1);
        fprintf(delete);
    end
else
    bstext = repmat('\b', 1, digits);
    tekst = [bstext '%0' num2str(digits) '.0f'];
    fprintf(tekst,i);
end




