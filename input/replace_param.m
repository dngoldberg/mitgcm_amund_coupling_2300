function replace_param(datafile,paramname,value,varargin)

if (length(varargin)==0)
    is_str=true;
else
    is_str=false;
end

if (~ischar(value));
    line1 = [' ' paramname ' = ' num2str(value) ','];
elseif (~is_str)
    line1 = [' ' paramname ' = ' value ','];
else
    line1 = [' ' paramname ' = ' strcat('''', value, '''') ','];
end
cmd = ['!sed "s/.*' paramname '.*/' line1 '/i" ' datafile ' > tmpfile; cp tmpfile ' datafile];
eval(cmd);

return