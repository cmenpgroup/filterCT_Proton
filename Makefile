.DELETE_ON_ERROR:

ifndef CLASTOOL
    $(error "Please set the variable CLASTOOL")
endif

ifndef ANALYSIS
    $(error "Please set the variable ANALYSIS")
endif

ROOTCONFIG := root-config
ROOTCFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTGLIBS := $(shell $(ROOTCONFIG) --glibs)

CXX := c++
CXXFLAGS := -O2 -Wall -fPIC $(ROOTCFLAGS)
LD := c++
LDFLAGS := -O2 $(ROOTLDFLAGS)

INCLUDES := -I$(ANALYSIS)/include \
               -I$(CLASTOOL)/include
LIBS := $(ROOTGLIBS) \
               -L$(CLASTOOL)/slib/${OS_NAME} -lClasTool \
               -L$(ANALYSIS)/slib/ -lTIdentificator

FILES := filterCT_Proton


.PHONY: all clean

all: $(FILES)

filterCT_Proton: filterCT_Proton.o


%: %.o
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@


clean:
	rm -f $(FILES:%=%.o) $(FILES)
