##
## @file		makefile
## @brief		MDFort  makefile.
## @author		Rinku Mishra <rinku.mishra@cppipr.res.in>
##

F95 = gfortran -w

LOCAL=local
INIFLAG=-I$(LOCAL)/include -L$(LOCAL)/lib -lcfgio
C_DIR := $(shell pwd)

# Definition of the Flags
OMPFLAGS = -O -fopenmp

EXEC	= mdfort

SDIR	= src
ODIR  = src/obj
LDIR	= lib
OUTDIR = output

SRC_ 	= main.f95# Additional CPP files
OBJ_	= $(SRC_:.f95=.o)

SRC = $(patsubst %,$(SDIR)/%,$(SRC_))
OBJ = $(patsubst %,$(ODIR)/%,$(OBJ_))

all: $(EXEC)

$(EXEC) : $(ODIR)/main.o $(OBJ)
	@echo "Linking MDFort"
	@$(F95) $^ -o $@ $(INIFLAG) $(OMPFLAGS)
	@echo "MDFort is built"
	@mkdir -p $(OUTDIR)

$(ODIR)/%.o: $(SDIR)/%.f95
	@echo "Compiling $<"
	@mkdir -p $(ODIR)
	@$(F95) -c $< -o $@ $(INIFLAG) $(OMPFLAGS)

subsystems:
	@cd $(LDIR)/iniparser && $(MAKE)
	export PATH=$(C_DIR)/$(LOCAL)/bin
	export C_INCLUDE_PATH=$(C_DIR)/$(LOCAL)/include
	export LIBRARY_PATH=$(C_DIR)/$(LOCAL)/lib

clean:
	@echo "Cleaning compiled files"
	@echo "run make veryclean to remove executables"
	@rm -f *~ $(ODIR)/*.o $(SDIR)/*.o $(SDIR)/*~
	@rm -rf $(OUTDIR)

veryclean: clean
	@echo "Cleaning executables and iniparser"
	@rm -f $(EXEC)
	@cd $(LDIR)/iniparser && $(MAKE) clean > /dev/null 2>&1

run:
	@echo "Running MDFort"
	./mdfort
