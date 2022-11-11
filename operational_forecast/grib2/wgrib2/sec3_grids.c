#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grb2.h"
#include "wgrib2.h"

/*
 * create the grib2 sec3 and latitudes and longitudes values for various grids
 *
 * public domain 6/2010 Wesley Ebisuzaki
 *
 */



#define ANGLE_FACTOR 1.0e6

/* create the grib2 sec3 and latitudes and longitudes for lat-lon grid */

unsigned char *sec3_lola(int nx, double x0, double dx, int ny, double y0, double dy, 
	unsigned char **old_sec, double **lon, double **lat) {

    static unsigned char gds[72];
    unsigned char *p;
    double r, *a, *b;
    int i, j;

    uint_char(72, gds);			/* 1-4 length of section */

    gds[4] = 3;				/* section 3 */
    gds[5] = 0;				/* source of grid  */

    uint_char(nx*ny, gds+6);		/* number of data points */
        
    gds[10] = 0;			/* number of octets for optional list of number */
    gds[11] = 0;			/* list */
    gds[12] = gds[13] = 0;		/* template 3.0  lat-lon grid*/

    /* shape of the earth */

    if (old_sec == NULL) {		/* use NCEP default */
	gds[14] = 6;			/* shape of earth */
	gds[15] = 0;			/* scale factor of radius of spherical earth */
	uint_char(0, gds+16);		/* scaled value of radius of spherical earth */
	gds[20] = 0;			/* scale factor of major radious of oblate earth */
	uint_char(0, gds+21);		/* scaled value of major axis */
	gds[25] = 0;			/* scale factor of minor axis */
	uint_char(0, gds+26);		/* scaled value of minor axis */
    }
    else {
	p = code_table_3_2_location(old_sec);
	for (i = 0; i < 16; i++) gds[14+i] = p[i];
    }

    uint_char(nx, gds+30);
    uint_char(ny, gds+34);
    uint_char(0, gds+38);		/* 0 */
    uint_char(0xffffffff, gds+42);	/* undefined */

    int_char((int) floor(y0*ANGLE_FACTOR+0.5), gds+46);        /* lat1 */
    uint_char((int) floor(x0*ANGLE_FACTOR+0.5), gds+50);        /* lon1 */

    gds[54] = 32+16;			/* flag table 3.3 */

    r = y0 + (ny-1)*dy;
    int_char((int) floor(r*ANGLE_FACTOR+0.5), gds+55);        /* lat1 */
    r = x0 + (nx-1)*dx;
    if (r >= 360.0) r -= 360.0;
    if (r < 0.0) r += 360.0;
    uint_char((int) floor(r*ANGLE_FACTOR+0.5), gds+59);        /* lon1 */

    int_char((int) floor(fabs(dx) * ANGLE_FACTOR+0.5), gds+63);        /* Di i direction increment */
    int_char((int) floor(fabs(dy) * ANGLE_FACTOR+0.5), gds+67);        /* Dj i direction increment */

    gds[71] = 0;                                /* scanning mode */
    if (dx < 0.0) gds[71] |= 128;
    if (dy > 0.0) gds[71] |= 64;

    /* now for the lat and lon values */
    if (lon == NULL || lat == NULL) return gds;
    *lon = (double *) malloc(nx * ny * sizeof(double));
    *lat = (double *) malloc(nx * ny * sizeof(double));
    if (*lon == NULL || *lat == NULL) fatal_error("sec3_lola: memory allocation","");

    a = *lat;
    b = *lon;;
    for (j = 0; j < ny; j++) {
	r  = j * dy + y0;
	for (i = 0; i < nx; i++) {
	    *b++ = i * dx + x0;
	    *a++ = r;
	}
    }
    return gds;
}

/*
 * create the grib2 sec3 for straight forward gaussian grid
 *
 * grib2 allows regional and thinned gaussian grids
 *  this routine only allows global gaussian grids
 */

unsigned char *sec3_gaussian(int nx, double x0, double dx, int ny, double y0, 
	unsigned char **old_sec, double **lon, double **lat) {

    static unsigned char gds[72];
    unsigned char *p, *tmp_sec[8];
    double *tlon, *tlat;
    double r;
    int i;

    if (x0 < 0) x0 += 360.0;
	
    uint_char(72, gds);                 /* 1-4 length of section */

    gds[4] = 3;                         /* section 3 */
    gds[5] = 0;                         /* source of grid  */

    uint_char(nx*ny, gds+6);            /* number of data points */

    gds[10] = 0;                        /* number of octets for optional list of number */
    gds[11] = 0;                        /* list */
    gds[12] = 0; gds[13] = 40;          /* template 3.0  gaussian grid*/

    /* shape of the earth */

    if (old_sec == NULL) {              /* use NCEP default */
        gds[14] = 6;                    /* shape of earth */
        gds[15] = 0;                    /* scale factor of radius of spherical earth */
        uint_char(0, gds+16);           /* scaled value of radius of spherical earth */
        gds[20] = 0;                    /* scale factor of major radious of oblate earth */
        uint_char(0, gds+21);           /* scaled value of major axis */
        gds[25] = 0;                    /* scale factor of minor axis */
        uint_char(0, gds+26);           /* scaled value of minor axis */
    }
    else {
        p = code_table_3_2_location(old_sec);
        for (i = 0; i < 16; i++) gds[14+i] = p[i];
    }

    uint_char(nx, gds+30);
    uint_char(ny, gds+34);
    uint_char(0, gds+38);               /* 0 */
    uint_char(0xffffffff, gds+42);      /* undefined */

    int_char((int) floor(y0*ANGLE_FACTOR+0.5), gds+46);        /* lat1 */
    uint_char((int) floor(x0*ANGLE_FACTOR+0.5), gds+50);        /* lon1 */

    gds[54] = 32+16;                    /* flag table 3.3 */

    r = -y0;
    int_char((int) floor(r*ANGLE_FACTOR+0.5), gds+55);        /* lat2 */
    r = x0 + (nx-1)*dx;
    if (r >= 360.0) r -= 360.0;
    if (r < 0.0) r += 360.0;
    uint_char((int) floor(r*ANGLE_FACTOR+0.5), gds+59);        /* lon2 */

    int_char((int) floor(fabs(dx) * ANGLE_FACTOR + 0.5), gds+63);        /* Di i direction increment */
    uint_char(ny/2, gds+67);

    gds[71] = 0;                                /* scanning mode */
    if (dx < 0.0) gds[71] |= 128;
    if (y0 < 0.0) gds[71] |= 64;

    if (lon == NULL || lat == NULL) return gds;

    tmp_sec[3] = gds;
    gauss2ll(tmp_sec, &tlat, &tlon);
    *lon = tlon;
    *lat = tlat;
    return gds;
}

/*
 * create the grib2 sec3 and latitudes and longitudes for mercator
 */

unsigned char *sec3_mercator(double lad, int nx, double x0, double dx, int ny, double y0, double dy,
        unsigned char **old_sec, double **lon, double **lat) {

    static unsigned char gds[100];
    unsigned char *p, *tmp_sec[8];
    double r;
    double *tlon, *tlat;
    int i;

    if (x0 < 0) x0 += 360.0;

    uint_char(72, gds);			/* 1-4 length of section */

    gds[4] = 3;				/* section 3 */
    gds[5] = 0;				/* source of grid  */
    uint_char(nx*ny, gds+6);		/* number of data points */
    gds[10] = 0;			/* number of octets for optional list of number */
    gds[11] = 0;			/* list */
    gds[12] = 0; gds[13] = 10;		/* template 3.0  mercator grid */

    /* shape of the earth */

    if (old_sec == NULL) {              /* use NCEP default */
        gds[14] = 6;                    /* shape of earth */
        gds[15] = 0;                    /* scale factor of radius of spherical earth */
        uint_char(0, gds+16);           /* scaled value of radius of spherical earth */
        gds[20] = 0;                    /* scale factor of major radious of oblate earth */
        uint_char(0, gds+21);           /* scaled value of major axis */
        gds[25] = 0;                    /* scale factor of minor axis */
        uint_char(0, gds+26);           /* scaled value of minor axis */
    }
    else {
        p = code_table_3_2_location(old_sec);
        for (i = 0; i < 16; i++) gds[14+i] = p[i];
    }
    uint_char(nx, gds+30);
    uint_char(ny, gds+34);
    int_char((int) floor(y0*ANGLE_FACTOR+0.5), gds+38);        /* lat1 */
    uint_char((int) floor(x0*ANGLE_FACTOR+0.5), gds+42);        /* lon1 */
    int_char((int) floor(lad*ANGLE_FACTOR+0.5), gds+47);
    r = y0 + (ny-1)*dy;
    int_char((int) floor(r*ANGLE_FACTOR+0.5), gds+51);        /* lat2 */
    r = x0 + (nx-1)*dx;
    if (r >= 360.0) r -= 360.0;
    if (r < 0.0) r += 360.0;
    uint_char((int) floor(r*ANGLE_FACTOR+0.5), gds+55);        /* lon2 */
    int_char((int) floor(0*ANGLE_FACTOR+0.5), gds+60);        /* only handle 0 orient */

    int_char((int) floor(fabs(dx) * ANGLE_FACTOR+0.5), gds+64);        /* Di i direction increment */
    int_char((int) floor(fabs(dy) * ANGLE_FACTOR +0.5), gds+68);        /* Dj i direction increment */

    gds[59] = 0;                                /* scanning mode */
    if (dx < 0.0) gds[59] |= 128;
    if (dy > 0.0) gds[59] |= 64;

    if (lon == NULL || lat == NULL) return gds;

    tmp_sec[3] = gds;
    mercator2ll(tmp_sec, &tlat, &tlon);
    *lon = tlon;
    *lat = tlat;
    return gds;
}

/*
 * create the grib2 sec3 and latitudes and longitudes for lambert conformal (conic) grid
 */

unsigned char *sec3_lc(double lov, double lad, double latin1, double latin2, int proj,
	int nx, double x0, double dx, int ny, double y0, double dy, 
	unsigned char **old_sec, double **lon, double **lat) {

    static unsigned char gds[81];
    int i;
    unsigned char *p, *tmp_sec[8];
    double *tlon, *tlat;

    if (x0 < 0) x0 += 360.0;

    uint_char(81, gds);			/* 1-4 length of section */

    gds[4] = 3;				/* section 3 */
    gds[5] = 0;				/* source of grid  */

    uint_char(nx*ny, gds+6);		/* number of data points */
        
    gds[10] = 0;			/* number of octets for optional list of number */
    gds[11] = 0;			/* list */

    gds[12] = 0;			/* gdt=30  lambert conformal */
    gds[13]= 30;

    /* shape of the earth */

    if (old_sec == NULL) {              /* use NCEP default */
        gds[14] = 6;                    /* shape of earth */
        gds[15] = 0;                    /* scale factor of radius of spherical earth */
        uint_char(0, gds+16);           /* scaled value of radius of spherical earth */
        gds[20] = 0;                    /* scale factor of major radious of oblate earth */
        uint_char(0, gds+21);           /* scaled value of major axis */
        gds[25] = 0;                    /* scale factor of minor axis */
        uint_char(0, gds+26);           /* scaled value of minor axis */
    }
    else {
        p = code_table_3_2_location(old_sec);
        for (i = 0; i < 16; i++) gds[14+i] = p[i];
    }

    uint_char(nx, gds+30);
    uint_char(ny, gds+34);

    int_char((int) floor(y0*ANGLE_FACTOR+0.5), gds+38);        /* lat1 */
    uint_char((int) floor(x0*ANGLE_FACTOR+0.5), gds+42);        /* lon1 */

    gds[46] = 32+16;						/* flag table 3.3 */
    int_char((int) floor(lad*ANGLE_FACTOR+0.5), gds+47);
    uint_char((int) floor(lov*ANGLE_FACTOR+0.5), gds+51);
    int_char( (int) floor(fabs(dx)*1000.+0.5), gds+55);        /* Di 1/1000 of a meter */
    int_char( (int) floor(fabs(dy)*1000.+0.5), gds+59);        /* Dj 1/1000 of a meter */
    gds[63] = proj;
    gds[64] = 0;                                /* scanning mode */
    if (dx < 0.0) gds[64] |= 128;
    if (dy > 0.0) gds[64] |= 64;
    int_char((int) floor(latin1*ANGLE_FACTOR+0.5), gds+65);
    int_char((int) floor(latin2*ANGLE_FACTOR+0.5), gds+69);
    int_char((int) (0*ANGLE_FACTOR), gds+73);		/* lat of southern pole ? */
    int_char((int) (0*ANGLE_FACTOR), gds+77);			/* lon of southern pole ? */

    if (lon == NULL || lat == NULL) return gds;

    tmp_sec[3] = gds;
    lambert2ll(tmp_sec, &tlat, &tlon);
    *lon = tlon;
    *lat = tlat;
    return gds;
}

/*
 * create the grib2 sec3 and latitudes and longitudes for polar stereographic
 */

unsigned char *sec3_polar_stereo(double lov, double lad, int proj,
        int nx, double x0, double dx, int ny, double y0, double dy,
        unsigned char **old_sec, double **lon, double **lat) {

    static unsigned char gds[65];
    int i;
    unsigned char *p, *tmp_sec[8];
    double *tlon, *tlat;

    if (x0 < 0) x0 += 360.0;

    uint_char(65, gds);                 /* 1-4 length of section */

    gds[4] = 3;                         /* section 3 */
    gds[5] = 0;                         /* source of grid  */

    uint_char(nx*ny, gds+6);            /* number of data points */

    gds[10] = 0;                        /* number of octets for optional list of number */
    gds[11] = 0;                        /* list */

    gds[12] = 0;                        /* gdt=20  polar stereographic */
    gds[13]= 20;

    /* shape of the earth */

    if (old_sec == NULL) {              /* use NCEP default */
        gds[14] = 6;                    /* shape of earth */
        gds[15] = 0;                    /* scale factor of radius of spherical earth */
        uint_char(0, gds+16);           /* scaled value of radius of spherical earth */
        gds[20] = 0;                    /* scale factor of major radious of oblate earth */
        uint_char(0, gds+21);           /* scaled value of major axis */
        gds[25] = 0;                    /* scale factor of minor axis */
        uint_char(0, gds+26);           /* scaled value of minor axis */
    }
    else {
        p = code_table_3_2_location(old_sec);
        for (i = 0; i < 16; i++) gds[14+i] = p[i];
    }

    uint_char(nx, gds+30);
    uint_char(ny, gds+34);
    int_char((int) floor(y0*ANGLE_FACTOR+0.5), gds+38);        /* lat1 */
    uint_char((int) floor(x0*ANGLE_FACTOR+0.5), gds+42);        /* lon1 */

    gds[46] = 32+16;			/* flag table 3.3 */

    int_char((int) floor(lad*ANGLE_FACTOR+0.5), gds+47);
    uint_char((int) floor(lov*ANGLE_FACTOR+0.5), gds+51);
    int_char((int) floor(fabs(dx)*1000+0.5), gds+55);        /* Di 1/1000 of a meter */
    int_char((int) floor(fabs(dy)*1000+0.5), gds+59);        /* Dj 1/1000 of a meter */
    gds[63] = proj;
    gds[64] = 0;                                /* scanning mode */
    if (dx < 0.0) gds[64] |= 128;
    if (dy > 0.0) gds[64] |= 64;
    if (lon == NULL || lat == NULL) return gds;

    tmp_sec[3] = gds;
    polar2ll(tmp_sec, &tlat, &tlon);
    *lon = tlon;
    *lat = tlat;
    return gds;
}
