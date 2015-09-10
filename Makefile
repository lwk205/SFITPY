# -*- makefile -*-
.PHONY: all clean

F2PY = f2py --fcompiler=gfortran
F2PY_FLAGS = --quiet

all: getpar.so getbmag.so getmag_wei.so

getpar.so: parcode/getpar.f
	${F2PY} ${F2PY_FLAGS} -m getpar -c $<
getbmag.so: bincode/getbmag.f
	${F2PY} ${F2PY_FLAGS} -m getbmag -c $<
getmag_wei.so: tricode/getmag_wei.f
	${F2PY} ${F2PY_FLAGS} -m getmag_wei -c $<

clean:
	${RM} *.so parcode/*.pyf bincode/*.pyf tricode/*.pyf

