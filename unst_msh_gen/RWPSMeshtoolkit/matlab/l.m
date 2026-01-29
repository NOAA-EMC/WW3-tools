function l(x)
if nargin<1
    ls -lrth
else
    eval(['ls ',x,' -lrth'])
end