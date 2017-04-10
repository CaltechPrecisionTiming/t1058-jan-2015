CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)

CPPFLAGS := $(shell root-config --cflags) -I$(INC)/include
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR) -lMathMore

CPPFLAGS += -g

TARGET  = SiPMvsPhotekNoAmp
TARGET1 = PlotTimeResolutionVsCharge

SRC  = app/SiPMvsPhotekNoAmp.cc src/tree.cc
SRC1 = app/PlotTimeResolutionVsCharge.cc

OBJ = $(SRC:.cc=.o)
OBJ1 = $(SRC1:.cc=.o)

all : $(TARGET) $(TARGET1)


$(TARGET) : $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET1) : $(OBJ1)
	$(LD) $(CPPFLAGS) -o $(TARGET1) $(OBJ1) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^


%.o : %.cc	
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@	
	@echo $<
clean :
	rm -f *.o src/*.o $(TARGET) $(TARGET1) *~

