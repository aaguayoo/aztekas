#####################################
#### define dimension parameters ####
#####################################
#----------------------------------------------------------------------
# From the list below, please activate/deactivate the options that     
# apply to your run. If you modify any of these options, make sure     
# that you recompile the whole code by typing "make clean; make" 
#----------------------------------------------------------------------
# one dimensional code
#PARAM += -Ddim=1
#PARAM += -Deq=3
#PARAM += -Dgraf=1
# two dimensional code
PARAM += -Ddim=2
PARAM += -Deq=4
PARAM += -Dgraf=2
# three dimensional code
#PARAM += -Ddim=3
#PARAM += -Deq=5
#PARAM += -Dgraf=3

# Cylindrical coordinates
PARAM += -Dalfa=1
# Cartesian coordinates
#PARAM += -Dalfa=1
# Spherical coordinates
#PARAM += -Dalfa=1
#####################################

SOURCES = main.c alloc.c array.c mesh.c initial.c \
          timestep.c input.c output.c integration.c flux.c \
          pvrs-rk.c limiters.c matrixes.c boundaries.c \
          ./HD/amatrix.c ./HD/dmmatrix.c ./HD/dnmatrix.c \
          ./HD/domatrix.c ./HD/qvector.c ./HD/fvector.c \
          ./HD/gvector.c ./HD/hvector.c ./HD/svector.c 

#OBJS = $(SOURCES:.c=.o)

FLAGS = -lm -Ofast

COMPILER = gcc

EXEC = aztekas

$(EXEC): $(SOURCES)
	@echo ""
	@echo "Compiling problem file ..."
	$(COMPILER) $(PARAM) -fopenmp $(SOURCES) $(FLAGS) -o $(EXEC)  
	@echo ""
	@echo "aztekas compiled successfully"
	
clean:
	rm -f $(EXEC)
	
