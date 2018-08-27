#####################################
#### define dimension parameters ####
#####################################
#----------------------------------------------------------------------
# From the list below, please activate/deactivate the options that     
# apply to your run. If you modify any of these options, make sure     
# that you recompile the whole code by typing "make clean; make" 
#----------------------------------------------------------------------

#####################################
# Tests
# 0. Custom problem
# 1. Riemann problem (1D, 2D, 3D)
# 2. Kelvin-Helmholtz (2D, 3D)
# 3. Jet (1D, 2D, 3D)
# 4. Spherical Accretion (1D, 2D, 3D)
# 5. Wind (1D, 2D, 3D)
#####################################

# Set test 
PARAM += -DTEST=4

#####################################
# Ghost cells
#####################################

PARAM += -Dgc=3

#####################################

#####################################
# Physics (rhd,hd,mhd,rmhd)
#####################################

PHY = hd

# -Grav if you want classical gravity
# nothing if you don't.
grav = -grav

#####################################

#####################################
# Dimension
#####################################

DIM = 2

ifeq ($(DIM),1)
	PARAM += -Ddim=1
	PARAM += -Deq=3
	PARAM += -Dgraf=1
endif
ifeq ($(DIM),2)
	PARAM += -Ddim=2
	PARAM += -Deq=4
	PARAM += -Dgraf=2
endif
ifeq ($(DIM),3)
	PARAM += -Ddim=3
	PARAM += -Deq=5
	PARAM += -Dgraf=3
endif

#####################################

#####################################
# Coordinates (cart,cyl,sph)
#####################################

COORD = sph

ifeq ($(COORD),cart)
	PARAM += -Dalfa=0
endif
ifeq ($(COORD),cyl)
	PARAM += -Dalfa=1
endif
ifeq ($(COORD),sph)
	PARAM += -Dalfa=2
endif

#####################################

#####################################
# Integration method
#####################################

#####################################

#####################################
# Compilation
#####################################

ifeq ($(PHY),hd)
SOURCES = main.c alloc.c array.c mesh.c initial.c initial_cond.c\
			 timestep.c input.c output.c integration.c flux.c \
			 rk-clas.c limiters.c vectors.c boundaries.c bound_cond.c custom.c\
			 ./HD-$(COORD)$(grav)/amatrix.c ./HD-$(COORD)$(grav)/dmmatrix.c \
			 ./HD-$(COORD)$(grav)/dnmatrix.c ./HD-$(COORD)$(grav)/domatrix.c  \
			 ./HD-$(COORD)$(grav)/q2uvector.c ./HD-$(COORD)$(grav)/u2qvector.c \
			 ./HD-$(COORD)$(grav)/qvector.c ./HD-$(COORD)$(grav)/fvector.c \
			 ./HD-$(COORD)$(grav)/gvector.c ./HD-$(COORD)$(grav)/hvector.c \
			 ./HD-$(COORD)$(grav)/svector.c ./HD-$(COORD)$(grav)/gauge.c 
endif

ifeq ($(PHY),rhd)
SOURCES = main.c alloc.c array.c mesh.c initial.c initial_cond.c\
			 timestep.c input.c output.c integration.c flux.c \
			 rk-rel.c limiters.c vectors.c boundaries.c bound_cond.c custom.c\
			 ./RHD-$(COORD)/amatrix.c ./RHD-$(COORD)/dmmatrix.c ./RHD-$(COORD)/dnmatrix.c \
			 ./RHD-$(COORD)/domatrix.c ./RHD-$(COORD)/q2uvector.c ./RHD-$(COORD)/u2qvector.c \
			 ./RHD-$(COORD)/qvector.c ./RHD-$(COORD)/fvector.c ./RHD-$(COORD)/gvector.c \
			 ./RHD-$(COORD)/hvector.c ./RHD-$(COORD)/svector.c ./RHD-$(COORD)/gauge.c  
endif

#OBJS = $(SOURCES:.c=.o)

FLAGS = -Ofast -lm

COMPILER = gcc

EXEC = aztekas

#####################################

$(EXEC): $(SOURCES)
	@echo ""
	@echo "Compiling problem file ..."
	$(COMPILER) $(PARAM) -fopenmp $(SOURCES) $(FLAGS) -o $(EXEC)  
	@echo "$(bla)"
	@echo "aztekas compiled successfully"
	
clean:
	rm -f $(EXEC)
	
