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




digits = ceil(log10(max(to,from) + 1));
% print first part if i is at the beguining or if it0 is out of interval
% and never when i is at the end
if (i == from ||  ~xor(it0<=from,it0<=to)) && (i ~= to)
    
    tekst = [label  num2str(max(to,from)) ':%0' num2str(digits) '.0f'];
    fprintf(tekst,i);
    % at the end delete everything if someting was printed before (it0 in the
    % range)
elseif i == to 
    if xor(it0<=from,it0<=to)
        delete = repmat('\b', 1, numel(label)+digits*2+1);
        fprintf(delete);
    end
    % remove last number and print the current
else
    bstext = repmat('\b', 1, digits);
    tekst = [bstext '%0' num2str(digits) '.0f'];
    fprintf(tekst,i);
end




