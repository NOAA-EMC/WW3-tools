%Make some plots of output Mesh and input data
SetPath
PlotJigsawGriddedData('Topo.DFun.NWcoastal.PSLG.msh',[-123.0885 -122.6996   48.3679   48.7568],[0,300]);
PlotJigsawUnstMesh('NWcoastal/RWPS.F.LLH.msh',[-123.0885 -122.6996   48.3679   48.7568],[0,3000],[0,8]);
PlotWW3Mesh('RWPS.WW3a1.msh',[-123.0885 -122.6996   48.3679   48.7568],[0,3000],[0,8]);
PlotJigsawPSLG('NWcoastal.PSLG.msh',[-123.0885 -122.6996   48.3679   48.7568]); % this can cause Out Of Memory errors for large scale PSLGs

%Below can cause OOM
%PlotJigsawGriddedData('DFun.NWcoastal.PSLG.msh',[-123.0885 -122.6996   48.3679   48.7568],[0,3*10^5]);
