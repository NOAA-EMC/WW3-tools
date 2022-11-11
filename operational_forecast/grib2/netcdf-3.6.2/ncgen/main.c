/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/ncgen/main.c,v 1.11 2006/08/29 15:26:04 russ Exp $
 *********************************************************************/

#include <config.h>
#include <stdio.h>		/* has getopt() under VMS */
#include <string.h>

#ifdef __hpux
#include <locale.h>
#endif
    
#include <netcdf.h>

#include "generic.h"
#include "ncgen.h"
#include "genlib.h"

extern int	yyparse(void);

const char *progname;			/* for error messages */
const char *cdlname;

int c_flag;
int fortran_flag;
int netcdf_flag;
int cmode_modifier;
int nofill_flag;
char *netcdf_name = NULL;	/* name of output netCDF file to write */

extern FILE *yyin;

static const char* ubasename ( const char* av0 );
static void usage ( void );
int main ( int argc, char** argv );


/* strip off leading path */
static const char *
ubasename(
	const char *av0)
{
	const char *logident ;
#ifdef VMS
#define SEP	']'
#endif
#ifdef MSDOS
#define SEP	'\\'
#endif
#ifndef SEP
#define SEP	'/'
#endif
	if ((logident = strrchr(av0, SEP)) == NULL)
		logident = av0 ;
	else
	    logident++ ;
	return logident ;
}


static void usage(void)
{
    derror("Usage: %s [ -b ] [ -c ] [ -f ] [ -k kind ] [ -x ] [ -o outfile]  [ file ... ]",
	   progname);
    derror("netcdf library version %s", nc_inq_libvers());
}


int
main(
	int argc,
	char *argv[])
{
    extern int optind;
    extern int opterr;
    extern char *optarg;
    int c;
    FILE *fp;

#ifdef __hpux
    setlocale(LC_CTYPE,"");
#endif
    
#ifdef MDEBUG
	malloc_debug(2) ;	/* helps find malloc/free errors on Sun */
#endif /* MDEBUG */

    opterr = 1;			/* print error message if bad option */
    progname = ubasename(argv[0]);
    cdlname = "-";

    c_flag = 0;
    fortran_flag = 0;
    netcdf_flag = 0;
    cmode_modifier = 0;
    nofill_flag = 0;

#if _CRAYMPP && 0
    /* initialize CRAY MPP parallel-I/O library */
    (void) par_io_init(32, 32);
#endif

    while ((c = getopt(argc, argv, "bcfk:l:no:v:x")) != EOF)
      switch(c) {
	case 'c':		/* for c output, old version of "-lc" */
	  c_flag = 1;
	  break;
	case 'f':		/* for fortran output, old version of "-lf" */
	  fortran_flag = 1;
	  break;
	case 'b':		/* for binary netcdf output, ".nc" extension */
	  netcdf_flag = 1;
	  break;
        case 'l':		/* specify language, instead of using -c or -f */
	    {
		char *lang_name = (char *) emalloc(strlen(optarg)+1);
		if (! lang_name) {
		    derror ("%s: out of memory", progname);
		    return(1);
		}
		(void)strcpy(lang_name, optarg);
		if (strcmp(lang_name, "c") == 0 || strcmp(lang_name, "C") == 0) {
		    c_flag = 1;
		}
		else if (strcmp(lang_name, "f77") == 0 || 
			 strcmp(lang_name, "fortran77") == 0 ||
			 strcmp(lang_name, "Fortran77") == 0) {
		    fortran_flag = 1;
		} else {	/* Fortran90, Java, C++, Perl, Python, Ruby, ... */
		    derror("%s: output language %s not implemented", 
			   progname, lang_name);
		    return(1);
		}
	    }
	  break;
	case 'n':		/* old version of -b, uses ".cdf" extension */
	  netcdf_flag = -1;
	  break;
	case 'o':		/* to explicitly specify output name */
	  netcdf_flag = 1;
	  netcdf_name = (char *) emalloc(strlen(optarg)+1);
	  if (! netcdf_name) {
	      derror ("%s: out of memory", progname);
	      return(1);
	  }
	  (void)strcpy(netcdf_name,optarg);
	  break;
	case 'x':		/* set nofill mode to speed up creation of large files */
	  nofill_flag = 1;
	  break;
        case 'v':		/* a deprecated alias for "kind" option */
	    /*FALLTHRU*/
        case 'k': /* for specifying variant of netCDF format to be generated */
	    {
		char *kind_name = (char *) emalloc(strlen(optarg)+1);
		if (! kind_name) {
		    derror ("%s: out of memory", progname);
		    return(1);
		}
		(void)strcpy(kind_name, optarg);
		/* The default kind is kind 1, with 32-bit offsets */
		if (strcmp(kind_name, "1") == 0 || 
		    strcmp(kind_name, "classic") == 0) {
		    cmode_modifier = 0;
		}
		/* The 64-bit offset kind (2)  should only be used if actually needed */
		else if (strcmp(kind_name, "2") == 0 || 
		    strcmp(kind_name, "64-bit-offset") == 0) {
		    cmode_modifier |= NC_64BIT_OFFSET;
		}
#ifdef USE_NETCDF4
		/* NetCDF-4 HDF5 format*/
		else if (strcmp(kind_name, "3") == 0 || 
		    strcmp(kind_name, "hdf5") == 0) {
		    cmode_modifier |= NC_NETCDF4;
		}
		/* NetCDF-4 HDF5 format, but using only nc3 data model */
		else if (strcmp(kind_name, "4") == 0 ||
		    strcmp(kind_name, "hdf5-nc3") == 0) {
		    cmode_modifier |= NC_NETCDF4 | NC_CLASSIC_MODEL;
		}
#endif 
		else 
		{
		   derror("Invalid format, try classic, 64-bit-offset, hdf5, or hdf5-nc3");
		   return 2;
		}
	    }
	  break;
	case '?':
	  usage();
	  return(8);
      }

    if (fortran_flag && c_flag) {
	derror("Only one of -c or -f may be specified");
	return(8);
      }

    argc -= optind;
    argv += optind;

    if (argc > 1) {
	derror ("%s: only one input file argument permitted",progname);
	return(6);
    }

    fp = stdin;
    if (argc > 0 && strcmp(argv[0], "-") != 0) {
	if ((fp = fopen(argv[0], "r")) == NULL) {
	    derror ("can't open file %s for reading: ", argv[0]);
	    perror("");
	    return(7);
	}
	cdlname = argv[0];
    }
    yyin = fp;
    return (yyparse());
}
