function geom=BoundaryShape2msh(S,flout)
% make coastlines for various smoothings of coastlines

%from example 6 aust.msh 

%geom = 
%    point: [1×1 struct]
%    edge2: [1×1 struct]
%    mshID: 'EUCLIDEAN-MESH'
%    fileV: 3

%geom.point.coord(1:10,:)
%  146.2929  -39.0150         0
%  146.2937  -39.0192         0
%  146.2846  -39.0242         0

%geom.edge2.index(1:10,:)
%           1           2           0
%           1       27577           0
%           2           3           0
           
%1km 
clear geom
geom.mshID='EUCLIDEAN-MESH'
geom.fileV = 3
if isstr(S)
    S = shaperead(S);
end

N=length(S);
x0=[];
y0=[];
edge0=[];
for k=1:N
    if mod(k,100)==0,k/N,end
    x=S(k).X(1:end-1);
    y=S(k).Y(1:end-1);
    if ~isempty(x)
        isisland=0;
        if and( x(1)==x(end),y(1)==y(end) )
            isisland=1;
        end
        if length(x)>1,
            if ~isisland,
                n0=length(x0);
                x0=[x0;x(:)];
                y0=[y0;y(:)];
                n1=length(x0);
                edge0=[edge0;[ n0+1:n1-1;n0+2:n1]'  ];
            else          
                n0=length(x0);
                x=x(1:end-1);
                y=y(1:end-1);          
                x0=[x0;x(:)];
                y0=[y0;y(:)];
                n1=length(x0);
                edge0=[edge0;[ n0+1:n1-1;n0+2:n1]';[n1,n0+1]];
            end
        end
    end
end

[ne,two]=size(edge0);
nn=length(x0);
geom.edge2.index=[edge0,zeros(ne,1)];
geom.point.coord=[x0(:),y0(:),zeros(nn,1)];

savemsh(flout,geom);
    