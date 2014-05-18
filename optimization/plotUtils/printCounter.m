function [] =  printCounter(from, to, i,label)
%PRINTCOUNTER diplaying counter in command view
%
% SYNTAX:
%   printCounter(FROM, TO, I) 
%   printCounter(To, I) %default from 1
%
%
if nargin < 3
    i = to;
    to = from;
    from = 1;
end
if nargin < 4
    label = [];
else
    label = [label ' ' ];
end

if to == -1
    digits = 8;
else
    digits = ceil(log10(to + 1));
end
if i == from
    if to == -1
        tekst = [label ':%0' num2str(digits) '.0f'];
    else
        tekst = [label  num2str(to) ':%0' num2str(digits) '.0f'];
    end
    fprintf(tekst,i);
elseif i == to
    delete = repmat('\b', 1, numel(label)+digits*2+1);
    fprintf(delete);
else
    bstext = repmat('\b', 1, digits);
    tekst = [bstext '%0' num2str(digits) '.0f'];
    fprintf(tekst,i);
end




