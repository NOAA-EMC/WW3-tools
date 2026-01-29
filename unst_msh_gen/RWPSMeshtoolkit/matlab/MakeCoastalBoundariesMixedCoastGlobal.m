function MakeCoastalBoundaries
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
clear 
close all
geom.mshID='EUCLIDEAN-MESH'
geom.fileV = 3

load GlobalCoastline1kmUSto15km.mat
N=length(S);
for k=1:N
    ns(k)=length(S(k).X(1:end-1));
end


close all ; initjig ;

% DEMO-2 -- generate a regionally-refined global grid with a
%   high-resolution "patch" (@37.5KM) embedded within a
%   uniform background grid (@150.KM).
%   JIGSAW is combined with a bisection procedure to improve
%   the regularity of the grid topology.

%------------------------------------ setup files for JIGSAW
rootpath='./'
  %  rootpath = fileparts( ...
  %      mfilename( 'fullpath' )) ;

    opts.geom_file = ...                % domain file
        fullfile(rootpath,...
            'cache','geom.msh') ;

    opts.jcfg_file = ...                % config file
        fullfile('./',...
            'cache','opts.jig') ;

    opts.mesh_file = ...                % output file
        fullfile(rootpath,...
            'cache','mesh.msh') ;

    opts.hfun_file = ...                % sizing file
        fullfile(rootpath,...
            'cache','spac.msh') ;

%------------------------------------ define JIGSAW geometry

    geom.mshID = 'ELLIPSOID-MESH' ;
    geom.radii = 6371 * ones(3,1) ;

   % opts.geom_file='geom'
    savemsh (opts.geom_file,geom) ;


%------------------------------------ compute HFUN over GEOM

%    topo = loadmsh(fullfile( ...
%        rootpath, 'files', 'topo.msh' ) );
    topoJS = loadmsh('/home/keston/jigsaw/jigsaw-geo-matlab/files/topo.msh');
    clear topo
    fl='RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc'
    lon=ncread(fl,'lon');
    lat=ncread(fl,'lat');
    z=ncread(fl,'bed_elevation');

       w=ones(3,3)/9;
       [nx,ny] = size(z)
       z3=conv2(z,w,'same');
       z3=z3(2:3:nx,2:3:ny);
      lon3=lon(2:3:nx);
      lat3=lat(2:3:ny);
    topo.point.coord{:,1}=lon3;
    topo.point.coord{:,2}=lat3;%fix to centers later
    topo.value=double(z3');


    topo.mshID='ELLIPSOID-GRID'
    topo.fileV=3;
     
    xpos = topo.point.coord{1};
    ypos = topo.point.coord{2};
    zlev = reshape( topo.value,length(ypos), length(xpos));

   [XPOS,YPOS] = meshgrid(xpos,ypos);

    XPOS = XPOS * pi/180 ;
    YPOS = YPOS * pi/180 ;

    S = shaperead('us_coastline/tl_2023_us_coastline.shp');
    xt=[];
    yt=[];
    for k=1:length(S);
        xt=[xt;S(k).X(1:end-1)'];
        yt=[yt;S(k).Y(1:end-1)'];
    end%1574701x1
    
 
    ZPOS=XPOS+i*YPOS;
     zt=xt+i*yt;
    [NX,NY]=size(ZPOS);
   D=zeros(NX,NY);
   t0=now;
    for k=1:NX,
        if mod(k,1)==0,k/NX,end
        for j=1:NY,
            z0=xpos(j)+i*ypos(k);
            d=min(abs(z0-zt));
             D(k,j)=d;
        end
        t1=now;
        disp(['eta:',num2str( 24*(NX-k) * (t1-t0) / k ),' hours'  ]);
    end
         
    save -v7.3 dist_D_HR3x3.mat D xpos ypos
     load dist_D_HR3x3.mat
    figure
    clf;pcolor(xpos,ypos,D);shading interp;
    colorbar;axis equal;colormap('jet');hold on
    plot(real(zt),imag(zt),'k.')

    figure
    hfun=+15. - 14.5* exp(-(D/2.5).^2);
    clf;pcolor(xpos,ypos,hfun);shading interp;
    colorbar;axis equal;colormap('jet');
    title('hfun')
    kprint('hfun.jpg')
    
    hmin=.5;%not nescesarry
    hfun(find(hfun<hmin))=hmin;
         
    save -v7.3 hfun3x3.mat hfun xpos ypos
  
    hmat.mshID = 'ELLIPSOID-GRID';
    hmat.radii = geom.radii;
    hmat.point.coord{1} = XPOS(1,:) ;
    hmat.point.coord{2} = YPOS(:,1) ;
    hmat.value = single(hfun);
    
   % opts.hfun_file='hfun'
    savemsh(opts.hfun_file,hmat) ;
   % opts.hfun_file='hfun.msh'
    hmat0 = marche(opts);%smooth hfun
 % no marche   hmat=hmat0;
%------------------------------------ build mesh via JIGSAW!

    fprintf(1,'  Constructing MESH...\n');

    opts.hfun_scal = 'absolute';
    opts.hfun_hmax = +inf ;
    opts.hfun_hmin = +0.0 ;

    opts.mesh_dims = +2 ;               % 2-dim. simplexes

    opts.optm_qlim = +.95 ;
    opts.optm_iter = + 32 ;

    rbar = mean(geom.radii(:)) ;
    hbar = mean(hmat.value(:)) ;

    nlev = round(log2( ...
        rbar / sin(.4 * pi) / hbar) ) ;

    mesh = tetris(opts, nlev - 1) ;

    plotsphere(geom,mesh,hmat,topo) ;

    drawnow ;
    set(figure(1),'units','normalized', ...
    'position',[.30,.50,.25,.30]) ;

    set(figure(2),'units','normalized', ...
    'position',[.05,.50,.25,.30]) ;
    set(figure(3),'units','normalized', ...
    'position',[.05,.10,.25,.30]) ;

    savemsh('mesh_mixed_coastline',mesh) ;

    meshP=GlobalMesh2Planer(mesh,topo);
    savemsh('mesh_mixed_coastline_planar',meshP) ;
    
    load GlobalCoastline1kmUSto15km.mat %land boiundary
    meshO = CutOutLand(meshP,S,3)
    savemsh('mesh_mixed_coastline_planar_ocean',meshO) ;

    
    save -v7.3 mesh_1kmUSCL_15kmOO_meshes.mat mesh meshP meshO
    meshOE = CutOutLandE(meshP,S)
    savemsh('mesh_mixed_coastline_planar_oceanE',meshOE) ;
