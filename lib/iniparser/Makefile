FC=gfortran
FFLAG=-O2
BIN=../../local/bin
INC=../../local/include
LIB=../../local/lib

$(shell mkdir -p $(BIN))
$(shell mkdir -p $(INC))
$(shell mkdir -p $(LIB))

TARGET=libcfgio.a
OBJECTS= string_conv_mod.o cfgio_mod.o

all: $(TARGET)

$(TARGET): $(OBJECTS)
	ar rc $@ $^
	ranlib $@
	cp *.mod $(INC)
	cp $(TARGET) $(LIB)

.SUFFIXES: .o .f90
%.o: %.f90
	$(FC) $(FFLAG) -c $< -I$(INC)

install: $(TARGET)
	# mkdir -p $(INC) $(LIB)
	cp *.mod $(INC)
	cp $(TARGET) $(LIB)

clean:
	rm *.mod *.o *.a
