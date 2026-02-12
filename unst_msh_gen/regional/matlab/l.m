function l(x)
%function l(str)
% just a shortcut for ls str -lrth from the matlab prompt
if nargin<1
    ls -lrth
else
    eval(['ls ',x,' -lrth'])
end
