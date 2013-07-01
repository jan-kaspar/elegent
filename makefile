GCC = g++ -Wall -O3 -fPIC
ROOTCINT = rootcint

# root includes and libs
INCS = `root-config --cflags` -I.
#LIBS = -L/usr/lib -lstdc++ `root-config --libs` -lMinuit
LIBS = `root-config --libs`

# my libraries source directories
MODDIR = IOMC/Elegent

# modules
modules = Model IslamModel PPPModel BSWModel BHModel JenkovszkyModel ExpModel InterpolationModel CoulombInterference Constants

#----------------------------------------------------------------------------------------------------

# directories to create
dirsToCreate = obj dict lib

modules_obj = $(addprefix obj/, $(addsuffix .o, $(modules)))
modules_dict_src = $(addprefix dict/, $(addsuffix _dict.cc, $(modules)))
modules_dict_obj = $(addprefix dict/, $(addsuffix _dict.o, $(modules)))

#----------------------------------------------------------------------------------------------------

.PHONY: all dirs dicts libs info clean

all: libs

libs: dirs dicts lib/libElegentModels.so

dirs: $(dirsToCreate)

$(dirsToCreate):
	mkdir -p $@

#dicts: $(modules_dict_obj)

$(modules_dict_obj) : dict/%_dict.o : dict/%_dict.cc
	$(GCC) $(INCS) -c $< -o $@

$(modules_dict_src) : dict/%_dict.cc : interface/%.h
	rm -f $@
	rm -f $(@:%.cc=%.h)
	$(ROOTCINT) $@ -c $<+

lib/libElegentModels.so: obj/Math.o $(modules_obj)
	@echo MAKING libElegentModels.so
	@$(GCC) -shared -olib/libElegentModels.so obj/*.o #dict/*.o

%.o : ../src/%.cc ../interface/%.h
	@echo COMPILING $<
	@$(GCC) $(INCS) -c $< -o $@

info:
	echo $(dirsToCreate)
	@echo
	echo $(modules)
	@echo
	echo $(modules_obj)
	@echo
	echo $(modules_dict_src)
	@echo
	echo $(modules_dict_obj)

clean:
	rm -rf obj
	rm -rf lib
	rm -rf dict
