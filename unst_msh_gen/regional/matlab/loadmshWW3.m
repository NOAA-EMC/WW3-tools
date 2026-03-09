function g=loadmshWW3(fl)
% function g=loadmshWW3(fl)
%
% load WW3 .msh format mesh see:
% https://polar.ncep.noaa.gov/waves/wavewatch/manual.v5.16.pdf
% 
% mesh structure g has fields:
%
%          x: longitute
%          y: latitude
%          z: bathymetric depth 
%          e: (ne x 3) element list
%      etype: element type
%      ntags: 
%    PhsEnty: 
%     EleNum: 
%    ElemUnk: 
%        bnd: open ocean boundary nodes
%    BndType: 
%    BndGeom: 

fp=fopen(fl,'r')

str=fgetl(fp);
while ~strcmp(str,'$Nodes'),
    str=fgetl(fp);
end
str=fgetl(fp);
nn=sscanf(str,'%i')

A=fscanf(fp,'%i %f %f %f\n',[4,nn])'; 
x=A(:,2);
y=A(:,3);
z=A(:,4);

%$EndNodes
%$Elements
%30559
%1  2  3  0  1  0  14443  14445  14360
%2  2  3  0  2  0  13218  13352  13350
%3  2  3  0  3  0  6821  6970  6969
%4  2  3  0  4  0  9927  10126  10084
%5  2  3  0  5  0  12886  12758  12885
%6  2  3  0  6  0  483  525  526
%7  2  3  0  7  0  9903  10095  10096
%8  2  3  0  8  0  5426  5301  5300
%9  2  3  0  9  0  13213  13346  13344

str=fgetl(fp);
str=fgetl(fp);
str=fgetl(fp);
ne=sscanf(str,'%i')
ke=0;
nbn=0;
bnd=[];
BndType=[];BndGeom=[];EleType=[];EleGeom=[];
A=fscanf(fp,'%i');
%assume 
if A(3)==2,
    j6i=A(3:6:length(A));
    j2=find(j6i==2);
    dj2=j2(2:end)-j2(1:end-1);
    jg=find(dj2==1);
    j=jg(end)+1;
    Abnd=A(1:j*6);
    Mbnd=reshape(Abnd,[6,j])';
    Aele=A(j*6+1:end);
    ne=length(Aele)/9;
    Mele=reshape(Aele,[9,ne])';
else %Assume no bnds
     Mbnd=[];
     ne=length(A)/9;
     Mele=reshape(A,[9,ne])';
end
if ~isempty(Mbnd)
    bnd=Mbnd(:,6);
    BndType=Mbnd(:,4);
    BndGeom=Mbnd(:,2);
else
    bnd=[];
    BndType=[];
    BndGeom=[];
end

ele=Mele(:,7:9);
etype=Mele(:,2);
ntags=Mele(:,3);
PhsEnty=Mele(:,4);
EleNum=Mele(:,5);
ElemUnk=Mele(:,6);


g.x=x;
g.y=y;
g.z=z;
g.e=ele;
g.etype=etype;
g.ntags=ntags; 
g.PhsEnty=PhsEnty;
g.EleNum=EleNum;
g.ElemUnk=ElemUnk;
g.bnd=bnd;
g.BndType=BndType;
g.BndGeom=BndGeom;

if(0)%read line by line
    for k=1:nn
        str=fgetl(fp);
        A=sscanf(str,'%i %f %f %f');
        node(k)=A(1);
        x(k)=A(2);
        y(k)=A(3);
        z(k)=A(4);
    end
    
    for k=1:ne
        if mod(k,ne)==100,
            disp(['reading element ',int2str(k),' of ',int2str(ne)]);
        end
        str=fgetl(fp);
        A=sscanf(str,'%i');
        if length(A)==9,
            ke=ke+1;
            etype(ke)=A(2);
            ntags(ke)=A(3);
            PhsEnty(ke)=A(4);
            EleNum(ke)=A(5);
            ElemUnk(ke)=A(6);
            ele(ke,:)=A([7,8,9])';
        end
        if length(A)==6,
            nbn=nbn+1;
            bnd(nbn)=A(6);
            BndType(nbn)=A(4);
            BndGeom(nbn)=A(2);
        end
    end
    g.x=x;
    g.y=y;
    g.z=z;
    g.e=ele;
    g.etype=etype;
    g.ntags=ntags;
    g.PhsEnty=PhsEnty;
    g.EleNum=EleNum;
    g.ElemUnk=ElemUnk;
    g.bnd=bnd;
    g.BndType=BndType;
    g.BndGeom=BndGeom;
 end

fclose(fp);


%Element connectivities
%Connectivities are described inside de Elements section, the section 
% surrounded by the$Elements and $EndElements tags. The First line contains
% the number of elements, the rest of lines contains several integers with 
% the following meaning:
%    1st integer: element numbering.
%    2nd integer: element type.
%    3rd integer: number of tags.
%    4th integer: physical entity; FEconv use this number to create PMH element groups.
%    5th integer: elementary entity; when saving a Gmsh mesh, FEconv also uses the element numbering as elementary entity.
%    Next integers: several node numbers that defines connectivities.
%MSH finite element types allowed in feconv
