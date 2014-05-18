function stringType = intType2stringType(type)
if ischar(type)
    stringType = type;
else
    if type == 1
        stringType = 'Injector';
    elseif type == -1
        stringType = 'Producer';
    else
        warning('check converter');
        stringType = type;
    end
end

end