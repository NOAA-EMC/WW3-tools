function WriteWW3MeshX(g,WW3FileOut)
% function WriteWW3Mesh(g,WW3FileOut)
% Output unstructured mesh to WW3 .msh format
%
% inputs:
%       g - finite element mesh structure with fields:
%           g.x (nn x 1) longitude coordinates of nodes
%           g.y (nn x 1) latitude coordinates of nodes
%           g.e (ne x 3) element matrix
%       WW3FileOut (output file name)
%

OpenBndNodes=FindOuterBndWW3(g,1)

x=g.x;y=g.y;z=g.z;e=g.e;

nn=length(x)
[ne,three]=size(e)
nb=length(OpenBndNodes);

fp=fopen(WW3FileOut,'w');
fprintf(fp,'$MeshFormat\n')
fprintf(fp,'2 0 8\n')
fprintf(fp,'$EndMeshFormat\n')
fprintf(fp,'$Nodes\n')
fprintf(fp,'%i\n',nn)

A=[[1:nn]' x(:),y(:),z(:)];

fprintf(fp,'%i %f %f %f\n',A')

fprintf(fp,'$EndNodes\n')
fprintf(fp,'$Elements\n')
fprintf(fp,'%i\n',ne+nb)

B=[[1:nb]', 15*ones(nb,1) ,2*ones(nb,1),zeros(nb,1),zeros(nb,1),OpenBndNodes(:)];

fprintf(fp,'%i %i %i %i %i %i\n',B')

C=[ nb+[1:ne]',2*ones(ne,1) ,3*ones(ne,1),zeros(ne,1),[1:ne]',zeros(ne,1),e];
fprintf(fp,'%i %i %i %i %i %i %i %i %i\n',C')
fprintf(fp,'$EndElements\n')
fclose(fp);



%76  2  3  0  1  0  77  76  1
%$MeshFormat
%2 0 8
%$EndMeshFormat
%$Nodes
%3070
%1  -72.0576782709  40.9902316949  4.2878041267
%2  -72.0521937363  40.9713426805  13.8250312805
%3  -72.0469687227  40.9523807323  20.9117679596

%3068  -72.585996  40.813074  1.747097373
%3069  -72.586352  40.818835  1.5
%3070  -72.589697  40.813418  1.5
%$EndNodes
%$Elements
%5855
%1  15  2  0  0  75
%2  15  2  0  0  74
%3  15  2  0  0  73

%73  15  2  0  0  3
%74  15  2  0  0  2
%75  15  2  0  0  1
%76  2  3  0  1  0  77  76  1
%77  2  3  0  2  0  76  2  1
%78  2  3  0  3  0  78  2  76

%5854  2  3  0  5779  0  3065  3069  3067
%5855  2  3  0  5780  0  3068  3067  3070
%$EndElements

