#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

extern int decode, flush_mode, latlon, file_append, save_translation;
extern double *lat, *lon;
extern int GDS_change_no;
extern int use_scale, dec_scale, bin_scale, wanted_bits, max_bits;
extern int npnts, nx, ny;
extern enum output_grib_type grib_type;

extern int msg_no;
/*
 * HEADER:100:irr_grid:output:3:make irregular grid, nearest neighbor, X=lon-lat list Y=radius (km) Z=output grib file
 */
int f_irr_grid(ARG3) {

    int i, k, m;
    struct local_struct {
        int ngrid;
        double *lon_lat_list, radius;
        FILE *out;
        int *iptr;
        int last_GDS_change_no;
    };
    struct local_struct *save;
    const char *t;
    unsigned char *p;
    double tmp;
    float *array;
    unsigned char *new_sec[8], *new_sec3;

    /* initialization phase */

    if (mode == -1) {
        decode = latlon = 1;
        *local = save = (struct local_struct *)malloc( sizeof(struct local_struct));
	if (save == NULL) fatal_error("irr_grid: memory allocation","");

	/* count number of colons */
	t = arg1;
	i = 0;
	while (*t) {
	    if (*t++ == ':') i++;
	}
	if (i % 2 != 1) fatal_error("irr_grid: need lon0:lat0:lon1:lat1:..:lonN:latN","");
        save->ngrid = (i + 1)/2;
	save->lon_lat_list = (double *) malloc(save->ngrid * 2 * sizeof(double));
	save->iptr = (int *) malloc(save->ngrid  * sizeof(int));
	if (save->lon_lat_list == NULL || save->iptr == NULL ) fatal_error("irr_grid: memory allocation","");

	t = arg1;
        k = sscanf(t, "%lf%n", &tmp, &m);
	if (k != 1) fatal_error("irr_grid: lat-lon list, %s",t);
	save->lon_lat_list[0] = tmp;
	t += m;
	for (i = 1; i < 2*save->ngrid; i++) {
            k = sscanf(t, ":%lf%n", &tmp, &m);
	    if (k != 1) fatal_error("irr_grid: lat-lon list, %s",t);
	    save->lon_lat_list[i] = tmp;
	    t += m;
	}

	for (i = 0 ; i < save->ngrid; i++) {
	    tmp = save->lon_lat_list[i*2];
	    if (tmp < 0.0) save->lon_lat_list[i*2] = tmp + 360.0;
	    if (tmp > 360.0) save->lon_lat_list[i*2] = tmp - 360.0;
	    if (fabs(save->lon_lat_list[i*2+1]) > 90.0) fatal_error("irr_grid: bad latitude","");
	}

	if (sscanf(arg2,"%lf",&(save->radius)) != 1) fatal_error("irr_grid: radius %s", arg2);
        if ((save->out = ffopen(arg3, file_append ? "ab" : "wb" )) == NULL) fatal_error("irr_grid could not open file %s", arg3);
    }
    if (mode < 0) return 0;

    save = (struct local_struct *) *local;

    if (save->last_GDS_change_no != GDS_change_no) {
        save->last_GDS_change_no = GDS_change_no;
        if (lat == NULL || lon == NULL || data == NULL) fatal_error("irr_grid: no val","");

        /* find the nearest points for the grid */
        closest_init(sec);
        for (i = 0; i < save->ngrid; i++) {
            save->iptr[i] = closest(sec, save->lon_lat_list[i*2+1], save->lon_lat_list[i*2]);
        }
    }

    array = (float *) malloc(save->ngrid * sizeof(float));
    new_sec3 = (unsigned char *) malloc((30+8*save->ngrid) * sizeof(unsigned char));
    if (array == NULL || new_sec3 == NULL) fatal_error("irr_grid: memory allocation","");

    /* sec3 = grid defintion */
    uint_char(30+save->ngrid*8, new_sec3);
    new_sec3[4] = 3;		// sec3
    new_sec3[5] = 0;		// use table 3.1
    uint_char(save->ngrid, new_sec3+6);
    new_sec3[10] = 0;		// no optional list octets
    new_sec3[11] = 0;
    uint2_char(130, new_sec3+12);

    p = code_table_3_2_location(sec);
    if (p == NULL) {  // no earth descripition
	for (i = 14; i < 30; i++) {
	    new_sec3[i] = 255;
	}
    }
    else {
	for (i = 14; i < 30; i++) {
	    new_sec3[i] = p[i-14];
	}
    }

    /* make new_sec[] with new grid definition */
    for (i = 0; i < 8; i++) new_sec[i] = sec[i];
    new_sec[3] = new_sec3;

    for (i = 0; i < save->ngrid; i++) {
	array[i] = save->iptr[i] >= 0 ?  data[save->iptr[i]] : UNDEFINED;
	int_char( (int) (save->lon_lat_list[i*2+1] * 1000000.0), new_sec3 + 30 + i*8);
	uint_char( (int) (save->lon_lat_list[i*2] * 1000000.0), new_sec3 + 34 + i*8);
    }

    grib_wrt(new_sec, array, save->ngrid, save->ngrid, 1, use_scale, dec_scale, bin_scale,
                wanted_bits, max_bits, grib_type, save->out);

    free(array);
    free(new_sec3);

    return 0;
}


/*
 * grid_130.c
 * Public Domain 6/2011 Wesley Ebisuzaki
 *
 * support for grid #130 
 *  irregular grid defined by lat-lon values
 *
 *  v1.0 6-2011
 */


int irr_grid2ll(unsigned char **sec, double **lat, double **lon) {
    unsigned char *gds;
    int i;
    double *llat, *llon;

    int nnx, nny, nres, nscan;
    unsigned int nnpnts;

    get_nxny(sec, &nnx, &nny, &nnpnts, &nres, &nscan);

    if (nnpnts == 0) return 0;

    if ((*lat = (double *) malloc(nnpnts * sizeof(double))) == NULL) {
        fatal_error("space_view2ll memory allocation failed","");
    }
    if ((*lon = (double *) malloc(nnpnts * sizeof(double))) == NULL) {
        fatal_error("space_view2ll memory allocation failed","");
    }
    llat = *lat;
    llon = *lon;
    gds = sec[3];

    for (i = 0; i < nnpnts; i++) {
	*llat++ = (double) 1e-6 * int4(gds+30+i*8);
	*llon++ = (double) 1e-6 * int4(gds+34+i*8);
    }
    return 0;
}

