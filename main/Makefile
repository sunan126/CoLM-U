# Makefile for CLM main program

include ../include/Makeoptions

DEF  = $(shell grep -i cama_flood ../include/define.h)
CaMa = $(word 1, ${DEF})

# Name of executable
EXECUTABLE = clmu.x

####################################################################
.DEFAULT :

OBJS =  precision.o                 \
	GlobalVars.o                \
	PhysicalConstants.o         \
	LC_Const.o                  \
	PFT_Const.o                 \
	ncio.o                      \
	MOD_PFTimeInvars.o          \
	MOD_PCTimeInvars.o          \
	MOD_UrbanTimeInvars.o       \
	MOD_TimeInvariants.o        \
	MOD_LuLccTimeInvars.o       \
	timemanager.o               \
	user_specified_forcing.o    \
	metdata.o                   \
	MOD_PFTimeVars.o            \
	MOD_PCTimeVars.o            \
	MOD_UrbanTimeVars.o         \
	MOD_TimeVariables.o         \
	MOD_LuLccTimeVars.o         \
	MOD_LuLccTMatrix.o          \
	GETMET.o                    \
	MOD_1D_PFTFluxes.o          \
	MOD_1D_PCFluxes.o           \
	MOD_1D_UrbanFluxes.o        \
	MOD_1D_Fluxes.o             \
	MOD_1D_Forcing.o            \
	MOD_2D_Fluxes.o             \
	MOD_2D_Forcing.o            \
	ThreeDCanopy.o              \
	ALBEDO.o                    \
	UrbanShortwave.o            \
	UrbanLongwave.o             \
	UrbanALBEDO.o               \
	Urbannetsolar.o             \
	ASSIM_STOMATA_conductance.o \
	FRICTION_VELOCITY.o         \
	LEAF_temperature.o          \
	LEAF_temperature_PC.o       \
	SOIL_thermal_parameters.o   \
	SOIL_SNOW_hydrology.o       \
	SNOW_Layers_CombineDivide.o \
	GLACIER.o                   \
	LAKE.o                      \
	SIMPLE_OCEAN.o              \
	CLMMAIN.o                   \
	eroot.o                     \
	flxwrite.o                  \
	groundfluxes.o              \
	groundtem.o                 \
	LAI_empirical.o             \
	LAI_readin.o                \
	LAI_readin_nc.o             \
	UrbanLAI_readin_nc.o        \
	LEAF_interception.o         \
	lpwrite.o                   \
	meltf.o                     \
	netsolar.o                  \
	rain_snow_temp.o            \
	newsnow.o                   \
	orb_coszen.o                \
	qsadv.o                     \
	monthly_in_situ_co2_mlo.o   \
	rd_forcing.o                \
	snowage.o                   \
	snowfraction.o              \
	THERMAL.o                   \
	tridia.o                    \
	vec2xy.o                    \
	wetbulb.o                   \
	UrbanFlux.o                 \
	UrbanGroundFlux.o           \
	UrbanRoofTem.o              \
	UrbanWallTem.o              \
	UrbanPerviousTem.o          \
	UrbanImperviousTem.o        \
	UrbanAnthropogenic.o        \
	UrbanHydrology.o            \
	UrbanTHERMAL.o              \
	UrbanIniTimeVar.o           \
	Urban_readin_nc.o           \
	IniTimeVar.o                \
	lakedepth_readin.o          \
	soil_parameters_readin.o    \
	HTOP_readin_nc.o            \
	LuLccInitialize.o           \
	LuLccEnergyConserve.o       \
	LuLccWaterConserve.o        \
	LuLccDRIVER.o               \
	UrbanCLMMAIN.o              \
	UrbanCLMDRIVER.o            \
	UrbanCLM.o

####################################################################

ifneq (${CaMa},\#define)
# Compile CoLM decoupled with river routing scheme (CaMa-Flood)

${EXECUTABLE} : ${OBJS}
	${FF} ${FOPTS} ${OBJS} -o $@ ${LDFLAGS}
	@echo 'making CLM completed!'

$(OBJS):%.o:%.F90 ../include/define.h
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

else
####################################################################

# The global river model CaMa-Flood (version 3.4.4)
CaMa_DIR = ../CaMa

# CaMa Flood Model modules directories
CaMa_MODS = -I$(CaMa_DIR)/mod \
            -I$(CaMa_DIR)/lib \
            -I$(CaMa_DIR)/src

# CaMa Flood Model libs (static) directories
CaMa_LIBS = $(CaMa_DIR)/src/srcs.a \
            $(CaMa_DIR)/mod/mods.a \
            $(CaMa_DIR)/lib/libs.a

#--------------------------------------------

${EXECUTABLE} : ${OBJS}
	${FF} ${FOPTS} ${OBJS} ${CaMa_LIBS} -o $@ ${LDFLAGS}
	@echo 'making CLM with CaMa Flood Model completed!'

$(OBJS):%.o:%.F90 ../include/define.h
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${CaMa_MODS}

endif
####################################################################

clean :
	${RM} -f *.o *.mod
