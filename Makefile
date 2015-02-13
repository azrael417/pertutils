# This file
MAKEFILE = Makefile

#----------------------------------------------------------------------
#  User choices - edit to suit 
#----------------------------------------------------------------------
# 0. Shell (if it really matters)

#SHELL = /bin/bash

#----------------------------------------------------------------------
# 1. Architecture

# Compiling for a parallel machine?  blank for a scalar machine
MPP = false
OMP = false
ARCH = MAC
INSTALLDIR=${HOME}/lib

#Compiler Type:
ifeq ($(ARCH),INTEL)
ifeq ($(strip ${OMP}),true)
  CC = icpc -openmp
else
  CC = icpc 
endif
endif

ifeq ($(ARCH),VAN)
ifeq ($(strip ${OMP}),true)
  CC = g++ -fopenmp
else
  CC = g++
endif
endif

ifeq ($(ARCH),CRAY)
CC = CC
endif

ifeq ($(ARCH),MAC)
ifeq ($(strip ${OMP}),true)
  CC = g++ -fopenmp
else
  CC = g++
endif
endif


#Compiler Optimization Level:
ifeq ($(ARCH),INTEL)
OPT              = -O2 -march=native -ipo
FLAGS            = -fPIC -Wall -debug -std=c++11
LD               = icpc
LINK             = -shared -Wl,-soname,libpertutils.so -o libpertutils.so.1.0 *.o
FINISH           = cp libpertutils.so.1.0 ${INSTALLDIR}/; ln -sf ${INSTALLDIR}/libpertutils.so.1.0 ${HOME}/lib/libpertutils.so
endif

ifneq (,$(filter $(ARCH),VAN CRAY))
OPT              = -O2 -march=native
FLAGS            = -fPIC -Wall -g
LD               = g++
LINK             = -shared -Wl,-soname,libpertutils.so -o libpertutils.so.1.0 *.o
FINISH           = cp libpertutils.so.1.0 ${INSTALLDIR}/; ln -sf ${INSTALLDIR}/libpertutils.so.1.0 ${INSTALLDIR}/libpertutils.so
endif

ifeq ($(ARCH),MAC)
OPT              = -O2 -march=native
FLAGS            = -fPIC -Wall -g
LD               = g++
LINK             = -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-compatibility_version,1.0,-current_version,1.0,-install_name,${INSTALLDIR}/libpertutils.dylib -o ${INSTALLDIR}/libpertutils.1.0.dylib *.o
FINISH           = ln -sf ${INSTALLDIR}/libpertutils.1.0.dylib ${INSTALLDIR}/libpertutils.dylib
endif


#Extra libraries and includes:
ifeq ($(ARCH),INTEL)
LIBADD = -lpthread -lm
MYINCLUDEDIR = -I./ -I../mathutils/
DEFINES = 
endif

ifneq (,$(filter $(ARCH),VAN CRAY))
LIBADD = -lpthread -lm
MYINCLUDEDIR = -I./ -I../mathutils/
DEFINES = 
endif

ifeq ($(ARCH),MAC)
LIBADD = -lpthread -lm
MYINCLUDEDIR = -I./ -I../mathutils/
DEFINES = 
endif

#Complete Set of Flags:
CFLAGS = ${DEFINES} ${MYINCLUDEDIR} ${FLAGS} ${OPT}

all::
	${CC} ${CFLAGS} -c gamma.cpp -o gamma.o ${LIBADD}
	${CC} ${CFLAGS} -c bartens.cpp -o bartens.o ${LIBADD}
	${CC} ${CFLAGS} -c runnings.cpp -o runnings.o ${LIBADD}
	${CC} ${CFLAGS} -c chiptfuncs.cpp -o chiptfuncs.o ${LIBADD}
	${LD} ${LINK}
	${FINISH}
clean::
	-/bin/rm -f *.o *.a *.so.*
