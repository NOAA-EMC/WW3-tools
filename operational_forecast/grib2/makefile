# REQUIRES GMAKE!!!!
#
# makefile for wgrib2
# 
# compiles every #@?! library needed by wgrib2
# then tries to compile wgrib2
#
# (1) must use gnu-make
# (2) the environment variable FC must be set to the fortran-90 or higher
#        in order to compile netcdf and the optional IPOLATES
# (3) the environment veriable CC is for the C compiler
#        often set by the system 
#
# mod 1/07 M. Schwarb (libgrib2c name change)
# mod 4/09 W. Ebisuzaki (use config.h)
# mod 6/10 W. Ebisuzaki ipolates
# mod 8/11 W. Ebisuzaki support environment variable FC={fortran 90+ compiler}
#              needed by optional netcdf4 and ipolates
#
#   Optional modules:
#
# NETCDF3: link in netcdf3 library to write netcdf3 files
#    change: USE_NETCDF3=1 and USE_NETCDF4=0 in configuration below
#
# NETCDF4: link in netcdf4 library to write netcdf3/4 files
#    change: USE_NETCDF3=0 and USE_NETCDF4=1 in configuration below
#    need to download netcdf4 and hdf5 libraries and to put into grib2 directory
#    need to define environment variable FC to be the command for the fortran compiler
#
# IPOLATES: link in IPOLATES library to interpolate to new grids
#    change: USE_IPOLATES=1 in configuration below
#    need to define environment variable FC to be the command for the fortran compiler
#    need to modify makefile and perhaps source code
#
#  MYSQL: link in interface to MySQL to write to mysql database
#    change: USE_MYSQL=1 in configuration below
#    need to have mysql installed
#    may need to modify makefile
#
#  UDF: add commands for user-defined functions and shell commands
#    change: USE_UDF=1 in configuration below
#
#  REGEX: use regular expression library, on by default
#    change: USE_REGEX=0 to turn off (configuration below)
#
#  TIGGE: ability for TIGGE-like variable names, on by default
#    change: USE_TIGGE=0 to turn off (configuration below)
#
#
# this version uses netcdf 3 libraries -- compile C only version
#  get from UCAR if compile doesn't work
#
#
#
# on NCEP AIX
# export CC=/usr/vacpp/bin/xlC_r
# export CPP=/usr/bin/cpp
# export FC=xlf
#
# for OS-X: uncomment line for makefile -f scripts/makefile.darwin
#
SHELL=/bin/sh

# 
# netcdf3: write netcdf files with netcdf-3 library
# netcdf4: write netcdf files with netcdf-4 library
# regex: regular expression package used by (match,not), POSIX-2
# tigge: enable -tigge option for tigge names
# mysql: write to mysql files
# ipolates: fortran interpolation library
# udf: user defined functions
#
# the flags are stored in wgrib2/config.h
#

# Warning do not set both USE_NETCDF3 and USE_NETCDF4 to one
USE_NETCDF3=1
USE_NETCDF4=0
USE_REGEX=1
USE_TIGGE=1
USE_MYSQL=0
USE_IPOLATES=0
USE_UDF=0

# often enviroment variable FC=fortran compiler, if not, define it here
# FC=gfortran

ifeq ($(USE_NETCDF3),1)
  ifeq ($(USE_NETCDF4),1)
    $(error ERROR, USE_NETCDF3 = 1 and USE_NETCDF4 = 1: can not link in 2 netcdf libraries)
  endif
endif

ifeq ($(USE_NETCDF4),1)
  ifeq ($(FC),)
    $(error ERROR, USE_NETCDF4 = 1: must set fortran90 compiler by environement variable FC)
  endif
endif

ifeq ($(USE_IPOLATES),1)
  ifeq ($(FC),)
    $(error ERROR, USE_IPOLATES = 1: must set fortran90 compiler by environement variable FC)
  endif
endif


# wCPPFLAGS has the directory of the includes 
# wLDFLAGS has the directory/name of the library

wCPPFLAGS:=-O3
wLDFLAGS:=
cwd:=${CURDIR}

CONFIG_H=${cwd}/wgrib2/config.h
a:=$(shell echo "/* config.h */" > ${CONFIG_H})

ifeq ($(USE_REGEX),1)
   a:=$(shell echo "\#define USE_REGEX" >> ${CONFIG_H})
else
   a:=$(shell echo "//\#define USE_REGEX" >> ${CONFIG_H})
endif

ifeq ($(USE_TIGGE),1)
   a:=$(shell echo "\#define USE_TIGGE" >> ${CONFIG_H})
else
   a:=$(shell echo "//\#define USE_TIGGE" >> ${CONFIG_H})
endif

# grib2c library

g:=${cwd}/g2clib-1.2.1
glib:=$g/libgrib2c.a
wLDFLAGS+=-L$g -lgrib2c
wCPPFLAGS+=-I$g

# Jasper

j=${cwd}/jasper-1.900.1
jlib=$j/src/libjasper/.libs/libjasper.a
wLDFLAGS+=-L$j/src/libjasper/.libs -ljasper
wCPPFLAGS+=-I$j/src/libjasper/include


ifeq ($(USE_NETCDF3),1)
   n:=${cwd}/netcdf-3.6.2
   nlib:=$n/libsrc/.libs/libnetcdf.a
   wLDFLAGS+=-L$n/libsrc/.libs -lnetcdf
   wCPPFLAGS+=-I$n/libsrc
   a:=$(shell echo "\#define USE_NETCDF3" >> ${CONFIG_H})
else
   a:=$(shell echo "//\#define USE_NETCDF3" >> ${CONFIG_H})
endif

ifeq ($(USE_NETCDF4),1)
   n4:=${cwd}/netcdf-4.1.1
   n4lib:=${n4}/libsrc4/.libs/libnetcdf.a
   h5:=${cwd}/hdf5-1.8.6
   h5lib:=${h5}/src/.libs/libhdf5.a
#   wLDFLAGS+=-L${n4}/libsrc/.libs -lnetcdf
   wLDFLAGS+=-L${n4}/libsrc4/.libs -lnetcdf -L${h5}/hl/src/.libs -lhdf5_hl -L${h5}/src/.libs -lhdf5
   wCPPFLAGS+=-I${n4}/libsrc4
   a:=$(shell echo "\#define USE_NETCDF4" >> ${CONFIG_H})
else
   a:=$(shell echo "//\#define USE_NETCDF4" >> ${CONFIG_H})
endif

ifeq ($(USE_MYSQL),1)
   wCPPFLAGS+=`mysql_config --cflags`
   wLDFLAGS+=`mysql_config --libs`
   a:=$(shell echo "\#define USE_MYSQL" >> ${CONFIG_H})
else
   a:=$(shell echo "//\#define USE_MYSQL" >> ${CONFIG_H})
endif

ifeq ($(USE_IPOLATES),1)
   ip:=${cwd}/iplib
   iplib:=${ip}/libipolate.a
   wLDFLAGS+=-L${ip} -lipolate 

#  for compiling with fortran library
#  wLDFLAGS+= (libraries need by the fortran code)
#  wCPPFLAGS+= -D(FORTRAN Name)   see New_grid.c

# for G95 - personal system
   ifeq ($(FC),g95)
#      wLDFLAGS+=-L/export/wesley/wd51we/g95-install/lib/gcc-lib/i686-unknown-linux-gnu/4.0.3 -lf95
      wLDFLAGS+=-L/export/wesley/wd51we/g95-install_64/lib/gcc-lib/x86_64-unknown-linux-gnu/4.0.3 -lf95
      wCPPFLAGS+=-DG95
    endif

# for gfortran - ubuntu and cygwin 1.7.7-1
   ifeq ($(FC),gfortran)
      wLDFLAGS+=-lgfortran
      wCPPFLAGS+=-DGFORTRAN
   endif

# for open64 fortran - personal system
   ifeq ($(FC),openf95)
      wLDFLAGS+=/export/wesley/wd51we/opt/x86_open64-4.2.5.1/lib/gcc-lib/x86_64-open64-linux/4.2.5.1/libfortran.a
      wLDFLAGS+=/export/wesley/wd51we/opt/x86_open64-4.2.5.1/lib/gcc-lib/x86_64-open64-linux/4.2.5.1/libffio.a
      wCPPFLAGS+=-DOPENF95
   endif

# NCEP CCS:
   ifeq ($(FC),xlf)
      wLDFLAGS+=-L/usr/lib - -lxlf90 
      wCPPFLAGS+=-DXLF
   endif

   a:=$(shell echo "\#define USE_IPOLATES" >> ${CONFIG_H})
else
   a:=$(shell echo "//\#define USE_IPOLATES" >> ${CONFIG_H})
endif

ifeq ($(USE_UDF),1)
   a:=$(shell echo "\#define USE_UDF" >> ${CONFIG_H})
else
   a:=$(shell echo "//\#define USE_UDF" >> ${CONFIG_H})
endif

# save fortran and C compiler names in config.h file

a:=$(shell echo "\#define CC \"${CC}\"" >> ${CONFIG_H})
a:=$(shell echo "\#define FORTRAN \"${FC}\"" >> ${CONFIG_H})

# png 

p=${cwd}/libpng-1.2.44
plib=$p/.libs/libpng.a
wLDFLAGS+=-L$p/.libs -lpng
wCPPFLAGS+=-I$p

# z

z=${cwd}/zlib-1.2.5
zlib=$z/libz.a
wLDFLAGS+=-L$z -lz
wCPPFLAGS+=-I$z

wLDFLAGS+=-lm
wCPPFLAGS+=-I/usr/include ${CPPFLAGS}

# -----------------------------------------------------

wLDFLAGS+=-lm
wCPPFLAGS+=-I/usr/include ${CPPFLAGS}

# -----------------------------------------------------

# check if make is GNU make else use gmake
make_is_gnu:=$(word 1,$(shell make -v))
ifeq ($(make_is_gnu),GNU)
   MAKE:=make
else
   MAKE:=gmake
endif




w=wgrib2
prog=$w/wgrib2

${prog}:        $w/*.c $w/*.h ${jlib} ${nlib} ${zlib} ${plib} ${h5lib} ${glib} ${n4lib} ${iplib}
	cd $w && export LDFLAGS="${wLDFLAGS}" && export CPPFLAGS="${wCPPFLAGS}" && ${MAKE}

fast:        $w/*.c $w/*.h ${jlib} ${nlib} ${zlib} ${plib} ${h5lib} ${glib} ${n4lib} ${iplib}
	cd $w && export LDFLAGS="${wLDFLAGS}" && export CPPFLAGS="${wCPPFLAGS}" && ${MAKE}


${jlib}:
	cp $j.tar.gz tmpj.tar.gz
	gunzip -n tmpj.tar.gz
	tar -xvf tmpj.tar
	rm tmpj.tar
	cd $j && ./configure --without-x --disable-libjpeg --disable-opengl && ${MAKE}

${plib}:	${zlib}
	cp $p.tar.gz tmpp.tar.gz
	gunzip -n tmpp.tar.gz
	tar -xvf tmpp.tar
	rm tmpp.tar
#       for OSX
#	export LDFLAGS="-L$z" && cd $p && export CPPFLAGS="${wCPPFLAGS}" && make -f scripts/makefile.darwin
#	for everybody else
	export LDFLAGS="-L$z" && cd $p && export CPPFLAGS="${wCPPFLAGS}" && ./configure --disable-shared && ${MAKE}

${zlib}:
	cp $z.tar.gz tmpz.tar.gz
	gunzip -n tmpz.tar.gz
	tar -xvf tmpz.tar
	rm tmpz.tar
	cd $z && ./configure && ${MAKE}


${glib}:	${jlib} ${plib} ${zlib}
	touch ${glib}
	rm ${glib}
	cd $g && export CPPFLAGS="${wCPPFLAGS}" && ${MAKE}


${nlib}:
	cp netcdf.tar.gz tmpn.tar.gz
	gunzip -n tmpn.tar.gz
	tar -xvf tmpn.tar
	rm tmpn.tar
	cd $n && ./configure --enable-c-only && ${MAKE} check

${n4lib}:	${zlib} ${h5lib}
	mkdir -p ${cwd}/zlib/include
	mkdir -p ${cwd}/zlib/lib
	cp ${z}/*.h ${cwd}/zlib/include/
	cp ${z}/*.a ${cwd}/zlib/lib/

	mkdir -p ${cwd}/hdf5/include
	cp  ${h5}/src/*.h ${h5}/hl/src/*.h ${cwd}/hdf5/include
	mkdir -p ${cwd}/hdf5/lib
	cp ${h5}/src/.libs/*a ${h5}/hl/src/.libs/*a ${cwd}/hdf5/lib/

	cp ${n4}.tar.gz tmpn.tar.gz
	gunzip -n tmpn.tar.gz
	tar -xvf tmpn.tar
	rm tmpn.tar
	cd ${n4} && ./configure --disable-fortran --disable-cxx --disable-dap --enable-netcdf-4 --with-zlib=${cwd}/zlib --with-hdf5=${cwd}/hdf5 && ${MAKE}

#	cd ${n4} && ./configure --enable-c-only && ${MAKE} check

#	cd ${n4} && ./configure --enable-netcdf-4 --enable-c-only --with-zlib=${cwd}/zlib --with-hdf5=${cwd}/hdf5 && ${MAKE} check

#	cd ${n4} && ./configure --enable-c-only && ${MAKE} check

#	cd ${n4} && ./configure --enable-netcdf-4 --enable-c-only && ${MAKE} check


${h5lib}:
	cp ${h5}.tar.gz tmph5.tar.gz
	gunzip -n tmph5.tar.gz
	tar -xvf tmph5.tar
	rm tmph5.tar
	cd ${h5} && ./configure --disable-shared --disable-fortran --with-zlib=$z && ${MAKE} && ${MAKE}
#	cd ${h5} && ./configure --disable-shared --disable-fortran --with-zlib=$z && ${MAKE} && ${MAKE} check

${iplib}:
	cd ${ip} && export F90=${F90} && make

clean:
	cd $w && ${MAKE} clean
	cd $g && touch junk.a junk.o && rm *.o *.a
	rm -rf $n
	rm -rf $j
	rm -rf $p
	rm -rf $z
	rm -rf ${n4}
	rm -rf ${h5}

