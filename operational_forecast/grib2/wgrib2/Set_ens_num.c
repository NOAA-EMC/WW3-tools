#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"
#include "CodeTable4_4.h"
/*
 * Set_ens_num.c
 *
 * converts PDT 0,1 -> 1,   8,11 -> 11
 *
 * 11/2011: Public Domain: Wesley Ebisuzaki
 *
 */


/*
 * HEADER:100:set_ens_num:misc:3:convert PDT 0,1 -> 1,  8,11 -> 11, X=code table 4.6 Y=pert num Z=num ens members
 */
int f_set_ens_num(ARG3) {

    int i, n, pdt;
    struct local_struct {
        int type_ens;
	int ens_fcst;
	int num_ens;
	int n_sec4;
	unsigned char *sec4;
    };
    struct local_struct *save;

    if (mode == -1) {
	*local = save = (struct local_struct *)malloc( sizeof(struct local_struct));
	i = atoi(arg1);
	if (i < 0 || i > 255) fatal_error("set_ens_num: code table 4.6 0..255 found %s", arg1);
	save->type_ens = i;
	i = atoi(arg2);
	if (i < 0 || i > 255) fatal_error("set_ens_num: pert_num 0..255 found %s", arg2);
	save->ens_fcst = i;
	i = atoi(arg3);
	if (i < 0 || i > 255) fatal_error("set_ens_num: num ens emembers 0..255 found %s", arg3);
	save->num_ens = i;

	save->n_sec4 = 0;
	save->sec4 = NULL;
    }
    save = *local;

    if (mode < 0) return 0;

    pdt = code_table_4_0(sec);


    if (pdt == 1 || pdt == 11) {
	sec[4][34] = save->type_ens;
	sec[4][35] = save->ens_fcst;
	sec[4][36] = save->num_ens;
	return 0;
    }
    if (pdt != 0 && pdt != 8) {
	fprintf(stderr,"set_ens_num: only works with product defn template 0,1,8 and 11\n");
	return 0;
   }

    n = GB2_Sec4_size(sec);
    if (n+3 > save->n_sec4) {
	save->n_sec4 = n + 3;
	if (save->sec4) free(save->sec4);
	save->sec4 = (unsigned char *) malloc(save->n_sec4 * sizeof(unsigned char) );
	if (save->sec4 == NULL) fatal_error("set_ens_num: memory allocation error","");
    }
    
    // now to add ensemble information
    for (i = 0; i < 34; i++) save->sec4[i] = sec[4][i];
    save->sec4[34] = save->type_ens;
    save->sec4[35] = save->ens_fcst;
    save->sec4[36] = save->num_ens;
    for (i = 34; i <= n; i++) save->sec4[i+3] = sec[4][i];

    uint_char(n+3, save->sec4);                   // length of sec[4]
    save->sec4[7] = 0;                            // pdt = 0 -> 1, 8 -> 11
    save->sec4[8] = pdt == 8 ? 11 : 1;
    sec[4] = save->sec4; 
    return 0;
}
