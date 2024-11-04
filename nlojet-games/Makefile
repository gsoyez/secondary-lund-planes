
# Makefile generated automatically by /home/gsoyez/work/jets/boosted-jets/common/scripts/mkcxx.pl ''
# run 'make make' to update it if you add new files

CXX = g++  # for macs - otherwise get c++ = clang
CXXFLAGS = -Wall -g -O2
INCLUDE =

HEPDIR=/home/gsoyez/work/HEP-software/install
NLOJETDIR=$(HEPDIR)/nlojet++-4.1.5

MKCXXDIR=/home/gsoyez/work/jets/boosted-jets/common/scripts

NLOJETLHAPDFDIR=$(HEPDIR)/lhapdf-6.5.3
NLOJETFASTJETDIR=$(HEPDIR)/fastjet
NLOJETCMDLINEDIR=/home/gsoyez/work/jets/boosted-jets/common/scripts/../helpers/CmdLine
NLOJETSIMPLEHISTDIR=/home/gsoyez/work/jets/boosted-jets/common/scripts/../helpers
NLOJETQCDDIR=/home/gsoyez/work/jets/boosted-jets/common/scripts/../analytic-bases

NLOJETINC=-O3,-Wall,-Wno-misleading-indentation,-I$(NLOJETFASTJETDIR)/include,-I$(NLOJETSIMPLEHISTDIR),-I$(NLOJETCMDLINEDIR),-I$(NLOJETLHAPDFDIR)/include,-I$(NLOJETQCDDIR),-I/home/gsoyez/work/jets/boosted-jets/common/scripts/../helpers,
NLOJETLD=-Wl,-L$(NLOJETLHAPDFDIR)/lib,-lLHAPDF -Wl,-L$(NLOJETFASTJETDIR)/lib,-lfastjet,-lfastjetcontribfragile -Wl,-L$(NLOJETCMDLINEDIR),-lCmdLine -Wl,-L$(NLOJETQCDDIR),-lanalytic-base,-lgsl,-lgslcblas
INCLUDE += $(LCLINCLUDE)

COMMONSRC = 
F77SRC = 
COMMONOBJ = 

PROGSRC = 
PROGOBJ = 

INCLUDE += 
LIBRARIES += 


all:   nlojet-3jets-born.nlojet 



make:
	/home/gsoyez/work/jets/boosted-jets/common/scripts/mkcxx.pl ''

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

realclean: clean
	rm -vf  

.cc.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.cpp.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.C.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.f.o:         $<
	$(F77) $(FFLAGS) -c $< -o $@

%.nlojet: %.cc 
	$(NLOJETDIR)/bin/create-nlojet-user $< $(COMMONSRC) -o $(subst .cc,,$<) -Wx,$(NLOJETINC) $(NLOJETLD)

depend:
	makedepend  $(LCLINCLUDE) -Y --   -- $(COMMONSRC) $(PROGSRC)
# DO NOT DELETE
