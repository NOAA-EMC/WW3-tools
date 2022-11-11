#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"


#ifdef USE_IPOLATES

/*
 * ncep_grids:
 *
 * part of the new_grid package
 *
 * this routine converts the calls with "ncep grid I" to a regular call
 *
 * understands t382, various ncep grid numbers
 */


void ncep_grids(const char **arg1, const char **arg2, const char **arg3) {

   int grid_no;
   if (strcmp(*arg1,"ncep") == 0 && strcmp(*arg2,"grid") == 0) {

	if (strcmp(*arg3,"t382") == 0) grid_no = 128;
	else if (strcmp(*arg3,"t62") == 0) grid_no = 98;
	else grid_no = atoi(*arg3);

	switch(grid_no) {
	    case 2:
		*arg1="latlon";
		*arg2="0:144:2.5";
		*arg3="90:73:-2.5";
		break;

	    case 3: 
		*arg1="latlon";
		*arg2="0:360:1";
		*arg3="90:181:-1";
		break;

	    case 98:
		*arg1 = "gaussian";
		*arg2 = "0:192:1.875";
		*arg3 = "88.542:94";
		break;

	    case 128:  /* t382 gaussian grid */
		*arg1 = "gaussian";
		*arg2 = "0:1152:0.3125";
		*arg3 = "89.760995:576";
		break;

	    case 173:
		*arg1 = "latlon";
		*arg2 = "0.041666666667:4320:0.083333333333";
		*arg3 = "89.95833333333:2160:-0.083333333333";
		break;

	    case 221:
		*arg1="lambert:253:50";
		*arg2="214.5:349:32463.41";
		*arg3="1:277:32463.41";
		break;

	    case 230:
		*arg1="latlon";
		*arg2="0:720:0.5";
		*arg3="90:361:-0.5";
		break;

	    case 249:
		*arg1="nps:210.0:60.0";
		*arg2="188.4:367:9867.89";
		*arg3="45.4:343:9867.89";
		break;
	}
    }
   return;
}
#endif
