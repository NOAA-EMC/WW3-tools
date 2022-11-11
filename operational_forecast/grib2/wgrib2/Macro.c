#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"


/* 3/2008 Public Domain Wesley Ebisuzaki
 * 3/2008 Manfred Schwarb added -V
 */

extern const char *item_deliminator;
extern int file_append, decode;

/*
 * HEADER:100:s:inv:0:simple inventory
 */

/*
 * this is a simple macro .. see how easy it is!
 * would be more complicated if functions used static variables
 * minor complication if need to set decode or latlon flags
 */

int f_s(ARG0) {

    if (mode >= 0) {
	f_t(CALL_ARG0);
	strcat(inv_out,item_deliminator);
	inv_out += strlen(inv_out);

	f_var(CALL_ARG0);
	strcat(inv_out,item_deliminator);
	inv_out += strlen(inv_out);

	f_lev(CALL_ARG0);
	strcat(inv_out,item_deliminator);
	inv_out += strlen(inv_out);

	f_ftime(CALL_ARG0);
	strcat(inv_out,item_deliminator);
	inv_out += strlen(inv_out);

	f_misc(CALL_ARG0);
    }
    return 0;
}

/*
 * HEADER:100:S:inv:0:simple inventory with minutes and seconds (subject to change)
 */

int f_S(ARG0) {

    if (mode >= 0) {
        f_T(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_var(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_lev(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_ftime(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_misc(CALL_ARG0);
    }
    return 0;
}


/*
 * HEADER:100:s_out:inv_output:1:simple inventory written to X
 */

int f_s_out(ARG1) {

    if (mode == -1) {
        if ((*local = (void *) ffopen(arg1,file_append ? "a" : "w")) == NULL)
                fatal_error("Could not open %s", arg1);
        return 0;
    }
    if (mode >= 0) {
	f_s(CALL_ARG0);
	fprintf((FILE *) *local, "%s", inv_out);
    }
    return 0;
}


/*
 * HEADER:100:verf:inv:0:simple inventory using verification time
 */
int f_verf(ARG0) {

    if (mode >= 0) {
        f_vt(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_var(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_lev(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_ftime(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

	f_misc(CALL_ARG0);
    }
    return 0;
}

/*
 * Manfred Schwarb
 */

/*
 * HEADER:100:V:inv:0:diagnostic output
 */

int f_V(ARG0) {
    int oldmode;
    if (mode == -1) decode = 1;
    if (mode >= 0) {
        f_vt(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_lev(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_ftime(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        oldmode=mode;
        mode=1;
        f_var(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);
        mode=oldmode;

        f_misc(CALL_ARG0);
        strcat(inv_out,"\n    ");
        inv_out += strlen(inv_out);

        f_stats(CALL_ARG0);
        strcat(inv_out,"\n    ");
        inv_out += strlen(inv_out);

        f_grid(CALL_ARG0);
        strcat(inv_out,"\n");
    }
    return 0;
}

