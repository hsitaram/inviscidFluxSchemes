SOURCEDIR = ../source
OBJFILES = ${SOURCEDIR}/constants.o ${SOURCEDIR}/inputdata.o ${SOURCEDIR}/limiters.o ${SOURCEDIR}/flux_utils.o ${SOURCEDIR}/fluxschemes.o\
	   ${SOURCEDIR}/ausm_schemes.o ${SOURCEDIR}/cusp_schemes.o   	${SOURCEDIR}/vars.o ${SOURCEDIR}/calc.o ${SOURCEDIR}/io.o

FORT=gfortran

#flags=-ffpe-trap=invalid,zero,overflow -g -static
FFLAGS = -O2 

%.o: %.f90
	${FORT} -c ${FFLAGS} $< -o $@

stube:	${OBJFILES}
	mv ${SOURCEDIR}/*.o .
	${FORT} ${FFLAGS} *.o ${SOURCEDIR}/main.f90 -o stube.exe

clean:
	rm -rf *.o ${SOURCEDIR}/*.o *.mod stube.exe
