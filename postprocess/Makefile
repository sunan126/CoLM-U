# Makefile for post processing

include ../include/Makeoptions

# name of executable
EXE1 = bin2netcdf

####################################################################
.DEFAULT :

OBJS =  precision.o \
	GlobalVars.o \
	PhysicalConstants.o \
	MOD_2D_Fluxes.o \
	bin2netcdf.o

####################################################################

${EXE1} : ${OBJS}
	${FF} ${FOPTS} ${OBJS} -o $@ ${LDFLAGS}
	@echo ''
	@echo 'Making bin2netcdf completed!'
	@echo ''

$(OBJS):%.o:%.F90 ../include/define.h
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

clean :
	${RM} -f *.o *.mod
