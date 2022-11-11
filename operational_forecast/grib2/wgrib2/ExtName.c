#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

/*
 * ext_name: extended variable names
 *
 * A one time, the variable name was sufficient to identify the field
 *  along came, probabilities (50% precent chance was different from a 10% chance)
 *              ensembles
 *              statistical processing
 *              mass_density and chemical type
 *
 * Now we have "compound" variables - ensembles of chemical-types
 * Sooner or later .. 30% chance, ensemble member, daily mean, O3
 *
 * To handle the current and future extentions
 *
 *  Part A:  -f_misc
 *           inventory to print out the extensions  
 *           format :A=value:B=value:C=value:
 *  Part B   getExtName
 *           like getName but returns extended name
 *
 * public domain 10/2010: Wesley Ebisuzaki
 */


/*
 * HEADER:100:misc:inv:0:variable name qualifiers like chemical, ensemble, probability, etc
 */
int f_misc(ARG0) {

    char *p;
    const char *string;
    char cbuffer[STRING_SIZE];
    int need_space = 0;
    int pdt, val;
    static int error_count = 0;

    if (mode < 0) return 0;

    pdt = GB2_ProdDefTemplateNo(sec);
    inv_out += strlen(inv_out);
    p = inv_out;

    *cbuffer = 0;
    inv_out = cbuffer;
    f_ens(CALL_ARG0);
    if (strlen(inv_out)) {
	if (need_space) strcat(p,":");
	strcat(p,inv_out);
	need_space = 1;
    }

    inv_out[0] = 0;
    f_prob(CALL_ARG0);
    if (strlen(inv_out)) {
	if (need_space) strcat(p,":");
	strcat(p,inv_out);
	need_space = 1;
    }

    inv_out[0] = 0;
    f_spatial_proc(CALL_ARG0);
    if (strlen(inv_out)) {
	if (need_space) strcat(p,":");
	strcat(p,inv_out);
	need_space = 1;
    }

    val = code_table_4_3(sec);
    if (val == 6 || val == 7) {
	if (need_space) strcat(p,":");
	strcat(p,"analysis/forecast error");
	need_space = 1;
    }

    if (pdt == 7) {
	if (need_space) strcat(p,":");
	strcat(p,"analysis/forecast error");
	need_space = 1;
    }
    else if (pdt == 10) {
	if (need_space) strcat(p,":");
	strcat(p,"percent probability");
	need_space = 1;
   }
   else if (pdt == 40) {
	if (need_space) strcat(p,":");
	strcat(p,"chemical=");
        val = code_table_4_230(sec);
        if (val >= 0) {
            if (GB2_MasterTable(sec) <= 4) {
                if (error_count++ <= 10) {
		    if (GB2_Center(sec) == 98) 
			fprintf(stderr, "Warning: possible incompatible chemistry table .. table turned off.\n");
		    else if (GB2_Center(sec) != 7)
			fprintf(stderr,
                    "Warning: if file made with ECMWF API, possible incompatible chemistry table\n");
		}
            }

            string = NULL;
            switch(val) {
#include "CodeTable_4.230.dat"
            }
            if (GB2_MasterTable(sec) <= 4 && GB2_Center(sec) == 98) {
                string = NULL;
            }
            if (string != NULL)  strcat(p,string);
            else {
		p += strlen(p);
		sprintf(p,"chemical_%d",val);
	    }
	}
    }
    return 0;
}

/*
 * HEADER:400:set_ext_name:misc:1:X=0/1 extended name on/off
 */

int use_ext_name = 1;

int f_set_ext_name(ARG1) {
    use_ext_name = atoi(arg1);
    return 0;
}

/*
 * HEADER:400:ext_name:inv:0:extended name, var+qualifiers
 */

int f_ext_name(ARG0) {
    if (mode >= 0) {
        getExtName(sec, mode, NULL, inv_out, NULL, NULL,".","_");
    }
    return 0;
}


/* 
  get extend name - need to change some characters
  space -> *space
  colon -> *space
  delim -> *delin
 */

int getExtName(unsigned char **sec, int mode, char *inv_out, char *name, char *desc, char *units,
   const char *delim, const char *space) {

    char *save_inv_out, misc_string[STRING_SIZE], *p;
    const char *s;
    int i, istart,iend, ieq;

    if (sec == NULL) return 1;

    /* arguments for ARG0 */
    float *data = NULL;
    unsigned int ndata = 0;
    void **local = NULL;

    getName(sec, mode, inv_out, name, desc, units);

    if (use_ext_name == 0) return 0;

    save_inv_out = inv_out;
    misc_string[0] = 0;
    inv_out = misc_string;
    f_misc(CALL_ARG0);
    inv_out = save_inv_out;

    if (name && misc_string != NULL && strlen(misc_string) != 0) {
	istart = 0;
        p = name + strlen(name);

        while (misc_string[istart] != 0) {
	    /* look for end of specification */
	    ieq = -1;
	    iend = istart + 1;
	    while (misc_string[iend] != 0 && misc_string[iend] != ':') {
		if (ieq == -1 && misc_string[iend] == '=') ieq = iend;
		iend++;
	    }

	    s = delim;
	    while (*s)  *p++ = *s++; 

	    if (ieq != -1) istart = ieq + 1;

	    for (i = istart; i < iend; i++) {
		if (misc_string[i] == ' ' || misc_string[i] == ':') {
		    s = space;
		    while (*s)  *p++ = *s++; 
		}
		else *p++ = misc_string[i];
	    }
	    istart = iend;
        }
        *p = 0;
    }
    return 0;
}

/*
 * HEADER:400:full_name:inv:0:extended name, var+qualifiers
 */

int f_full_name(ARG0) {
    if (mode >= 0) {
        getExtName(sec, mode, NULL, inv_out, NULL, NULL,".","_");
	inv_out += strlen(inv_out);
	*inv_out++ = '.';
	*inv_out = 0;
	f_lev0(CALL_ARG0);
    }
    return 0;
}
