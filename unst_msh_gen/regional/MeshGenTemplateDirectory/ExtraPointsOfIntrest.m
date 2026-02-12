function [xus,yus]=ExtraPointsOfIntrest
xus=[];
yus=[];

%Now add Pacific Territories and COFA points
%------------------------------------------------------------------------
%NWS Pacific Region, via the Compact of Free Association, 
% %oversees operations of 5 Weather Service Offices (WSO) across 
% Micronesia in the Republic of Palau, Federated States of Micronesia 
% and the Republic of the Marshall Islands. WFO Guam provides routine
% forecasts as well as WWA services for these areas. -Eric Lau
%Palau,Yap, Chuuk, Pohnpei, Majuro, Pago Pago, Wake island

%Wake Island
coord=[19.2796,166.6499];%E
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Palau
coord=[7.4942, 134.5690]
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Yap:9.5557° N, 138.1399° E
coord=[9.5557, 138.1399]
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Chuuk: 7°25′N 151°47′E
coord=[7.374227, 151.754606]%E
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Pohnpei › Coordinates: 6.8519° N, 158.2147° E
coord=[6.8519, 158.2147]
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Majuro › Coordinates7.0667° N, 171.2667° E
coord=[7.0667, 171.2667]
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Pago Pago › Coordinates: 14.2732° S, 170.7030° W
coord=[-14.2732, -170.7030]% SW
xus=[xus,coord(2)];yus=[yus,coord(1)];
% Kosrae › Coordinates :5.3096° N, 162.9815° E
coord=[5.3096, 162.9815]
xus=[xus,coord(2)];yus=[yus,coord(1)];
%------------------------------------------------------------------------

%Marshall Islands
%Majuro
coord=[7.0667, 171.2667]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
%Ebeye
coord=[8.7815, 167.7373]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
%Micronesia
%Kolonia,
coord=[6.9636, 158.2102]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
%Pohnpei,
coord=[6.8519, 158.2147]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
% Chuuk-Weno
coord=[7.4523, 151.8422]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
% Tofol
coord=[5.3256, 163.0086]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
%Colonia -between Palau and Guam
coord=[9.5164,138.1222]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
n1=length(xus);
%------------------------------------------------------------------------
%Points from Curt
%Howland Island -Baker
coord=[0.8113, -176.6183]%W
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Johnston atoll
coord=[16.7295, -169.5336]%W
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Palmyra
coord=[5.8885, -162.0787]%W
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Jarvis Island
coord=[0.3744, -159.9967]%W
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Baker Island
%0.1936° N, 176.4769° W
coord=[0.1936,-176.4769 ]%W
xus=[xus,coord(2)];yus=[yus,coord(1)];
%------------------------------------------------------------------------
 coord=[18.4101, -75.0115 ]%W
 xus=[xus,coord(2)];yus=[yus,coord(1)];
 %If your team is going to the effort to add resolution for Majuro in RMI, I would strongly suggest doing the same for these atolls in RMI:
 %1. Kwajalein Atoll, which is home to part of the Ronald Regard Ballistic Missile Test Site (https://home.army.mil/kwajalein/index.php) and the underprediction of wave heights by NWS' WaveWatchIII model in January 2024, which I discussed on our last call, is what led to significant damage at the base, per: https://www.youtube.com/shorts/jH-pGoQDdcg   
 %2. Enewetak Atoll, which is the home of the Runit Dome (https://en.wikipedia.org/wiki/Runit_Island), is threatened by wave-driven overwash that has serious implications for the US Department of State, Department of Defense/War, and Department of the Interior via the Intergovernmental Compact of Free Association (https://www.doi.gov/oia/compacts-of-free-association).
 %3. Bikini Atoll, for similar reasons as Enewetak, although the radionuclides are all over the place and not all dumped in one location.
 %I would note that most of the atolls drop off at a 70-80 degree slope from approximately 30 m depth (which is generally less than 1000 m from shore) to over 1000 m depth, so there is not a need for a large region of increasing resolution to capture a broad continental shelf as characterizes CONUS.
 
 %Kwajalein Atoll   - 9.1898° N, 167.4243° E
 %Enewetak Atoll     - 11.4654° N, 162.1890° E
 %Bikini Atoll       - 11.6065° N, 165.3768° E
 
 %Kwajalein Atoll   - 9.1898° N, 167.4243° E
 coord=[9.1898, 167.4243 ]%W
 xus=[xus,coord(2)];yus=[yus,coord(1)];
 %Enewetak Atoll     - 11.4654° N, 162.1890° E
 coord=[ 11.4654, 162.1890]%W
 xus=[xus,coord(2)];yus=[yus,coord(1)];
 %Bikini Atoll       - 11.6065° N, 165.3768° E
 coord=[11.6065, 165.3768 ]%W
 xus=[xus,coord(2)];yus=[yus,coord(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add points to off shore banks we want to refine here we have Georges Bank
%and banks around the Bahamas.

xFB= [ -79.9909  -78.1963  -78.3957]
yFB=  [23.7978   26.8928   24.1530]
xGB = -67.4517
yGB =   41.3567

xus=[xus(:);xFB(:);xGB(:)]';
yus=[yus(:);yFB(:);yGB(:)]';

xus=LonCon(xus);
