# Makefile for the ROOT programs

ROOTCONFIG   := root-config
IncDir        = include/
SrcDir        = src/
BinDir        = bin/
LibDir        = lib/
HdrSuf        = h
ObjSuf        = o
SrcSuf        = cc
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

ARCH         := $(shell $(ROOTCONFIG) --arch)
PLATFORM     := $(shell $(ROOTCONFIG) --platform)

ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTINCDIR   := $(shell $(ROOTCONFIG) --incdir)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
HASTHREAD    := $(shell $(ROOTCONFIG) --has-thread)
ROOTDICTTYPE := $(shell $(ROOTCONFIG) --dicttype)
ROOTCINT     := rootcint
LIBNAME	      = SiPM
OPT2          = -g
CXX          := $(shell $(ROOTCONFIG) --cxx)
LD           := $(shell $(ROOTCONFIG) --ld)
CXXFLAGS      = $(OPT2) -Wall -fPIC
CXXFLAGS     += -I$(IncDir) $(ROOTCFLAGS)
SOFLAGS       = -shared
#MYLDFLAGS     = -L. -l$(LIBNAME)
FFTWLIBDIR    = /usr/local/lib
EXTRAROOTLIBS = -lSpectrum -lFFTW
LDFLAGS       =  $(OPT2)
LDFLAGS       += $(ROOTLDFLAGS)
LDFLAGS       += $(MYLDFLAGS)
LDFLAGS       += $(ROOTLIBS)
LDFLAGS       += $(EXTRAROOTLIBS)
SOLDFLAGS       =  $(OPT2)
SOLDFLAGS       += $(ROOTLDFLAGS)
SOLDFLAGS       += $(ROOTLIBS)
SOLDFLAGS       += $(EXTRAROOTLIBS)


EMPTY   := 
CLASSES_H  := $(shell ls $(IncDir)*.$(HdrSuf) | grep -v Dict)
CLASSES_S  := $(subst $(IncDir),$(SrcDir),$(CLASSES_H))
CLASSES_S  := $(subst .$(HdrSuf),.$(SrcSuf),$(CLASSES_S))
CLASSES_O  := $(subst .$(SrcSuf),.$(ObjSuf),$(CLASSES_S))
CLASSES_LDH:= $(subst .$(HdrSuf),LinkDef.$(HdrSuf),$(CLASSES_H)) #LINK DEF
CLASSES_DH := $(subst .$(HdrSuf),Dict.$(HdrSuf),$(CLASSES_H))
CLASSES_DS := $(subst .$(SrcSuf),Dict.$(SrcSuf),$(CLASSES_S))
CLASSES_DO := $(subst .$(SrcSuf),Dict.$(ObjSuf),$(CLASSES_S))
MYLIB   = $(LibDir)lib$(LIBNAME).$(DllSuf)

SRC     := $(shell ls $(SrcDir)ana*.$(SrcSuf))
SRC     += $(shell ls $(SrcDir)macro*.$(SrcSuf))
SRC     += $(shell ls $(SrcDir)get*.$(SrcSuf))
#SRC     += gainana.C getgain.C getdcr.C
OBJS 	:= $(subst .$(SrcSuf),.$(ObjSuf),$(SRC))
PROGS	:= $(subst .$(SrcSuf),,$(SRC))
PROGS	:= $(subst $(SrcDir),$(BinDir),$(PROGS))

.PHONY: all clean install lib

.SUFFIXES: .o .C

.C.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

all: $(PROGS) #$(MYLIB)

lib: $(MYLIB)

install: $(MYLIB)
	cp $(MYLIB) /usr/lib

##COMMENTED OUT BECAUSE IT USES THE IMPLICIT RULE FOR SINGLE OBJECT FILES
$(BinDir)%: $(SrcDir)%.o $(CLASSES_O)# $(CLASSES_DO) 
	$(LD) $(CXXFLAGS) $^ $(LDFLAGS) $(OutPutOpt) $@

clean:
	echo cleaning objs, dicts, progs and lib
	rm -fr $(OBJS) $(CLASSES_O) $(CLASSES_DO) $(CLASSES_DH) $(CLASSES_LDH) $(CLASSES_DS) $(PROGS) $(MYLIB)

%Dict.C %Dict.h: %.C %.h %LinkDef.h
	@echo "Generating dictionary for $*"
	rootcint -f $*Dict.C -c $*.h $*LinkDef.h  

%LinkDef.h : %.h
	@rm -f $@
	@echo "Generating $@"
	@echo "#ifdef __CINT__" > $@
	@echo "" >> $@
	@echo "#pragma link off all globals;" >> $@
	@echo "#pragma link off all classes;" >> $@
	@echo "#pragma link off all functions;" >> $@
	@echo "" >> $@
	@echo "#pragma link C++ class $*;" >> $@
	@echo "" >> $@
	@echo "#pragma link C++ global gROOT;" >> $@
	@echo "" >> $@
	@echo "#endif" >> $@

$(MYLIB): $(CLASSES_O) # $(CLASSES_DO)
	$(LD) $(SOFLAGS) $(SOLDFLAGS) $^ $(OutPutOpt) $@
