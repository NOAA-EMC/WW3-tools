#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "wgrib2.h"

/*
 * a simple extension to fopen
 *
 * if file is already open, just return the pointer to the file handler
 *
 * public domain 2/2008 Wesley Ebisuzaki
 *  v1.1 WNE: upgrade: add read/write detection
 *		- can be stdin and stdout (in previous version)
 *		if pipe - set flush_mode, error checking
 *          allow - to be stdin
 *
 */

extern int flush_mode;


FILE *ffopen(const char *filename, const char *mode) {

    struct opened_file { char *name; FILE *handle; int is_read_file; 
	struct opened_file *next; };

    static struct opened_file *opened_file_start = NULL;
    struct opened_file *ptr;
    struct stat stat_buf;  /* test for pipes */
    int is_read_file;
    const char *p;

    /* see if is a read/write file */
    is_read_file = 0;
    p = mode;
    while (*p) {
	if (*p++ == 'r') is_read_file = 1;
    }

    if (strcmp(filename,"-") == 0) {
	if (is_read_file) return stdin;
	flush_mode = 1;
//fprintf(stderr,"\n flush mode set by -\n");
	return stdout;
    }

    ptr = opened_file_start;
    while (ptr != NULL) {
	if (strcmp(filename,ptr->name) == 0) {
	    if (is_read_file != ptr->is_read_file) fatal_error("ffopen: file can only read or write not both: %s",
		ptr->name);
	    return ptr->handle;
	}
	ptr = ptr->next;
    }

    ptr = (struct opened_file *) malloc( sizeof(struct opened_file) );
    ptr->handle = fopen(filename,mode);
    if (ptr->handle == NULL) fatal_error("ffopen: could not open %s", filename);
    ptr->name = (char *) malloc(strlen(filename) + 1);
    strcpy(ptr->name, filename);
    ptr->is_read_file = is_read_file;
    ptr->next = opened_file_start;
    opened_file_start = ptr;

    /* check if output is to a pipe */
    if (is_read_file == 0) {
        if (stat(filename, &stat_buf) == -1) fatal_error("ffopen: could not stat file: %s", filename);
        if (S_ISFIFO(stat_buf.st_mode)) {
	    flush_mode  = 1;
//fprintf(stderr,"\n flush mode set by %s\n", filename);
	}
    }

    return ptr->handle;
}
