##
## @file		makefile
## @brief		MDFort  makefile.
## @author		Rinku Mishra <rinku.mishra@cppipr.res.in>
##

F95 = ~/opt/miniconda3/bin/gfortran

# Definition of the Flags
OMP = -O -fopenmp
#FFTWI = -I/usr/bin/include
#FFTWL = -lfftw3 -lm

EXEC	= mdfort

SDIR	= src
ODIR  = src/obj
OUTDIR = output

SRC_ 	= # Additional CPP files
OBJ_	= $(SRC_:.f95=.o)

SRC = $(patsubst %,$(SDIR)/%,$(SRC_))
OBJ = $(patsubst %,$(ODIR)/%,$(OBJ_))

all: $(EXEC)

$(EXEC) : $(ODIR)/main.o $(OBJ)
	@echo "Linking MDFort"
	@$(F95) $^ -o $@ $(OMP)
	@echo "MDFort is built"
	@mkdir -p $(OUTDIR)

$(ODIR)/%.o: $(SDIR)/%.f95
	@echo "Compiling $<"
	@mkdir -p $(ODIR)
	@$(F95) -c $< -o $@ $(OMP)

clean:
	@echo "Cleaning compiled files"
	@rm -f *~ $(ODIR)/*.o $(SDIR)/*.o $(SDIR)/*~
	@rm -rf $(OUTDIR)

# data:
# 	@echo "moving data files"
# 	@mv fort.* $(OUTDIR)/
