# ROOT includes and libs
INCS = `root-config --cflags` -I.
LIBS = `root-config --libs`

# HepMC includes and libs
INCS += -I/usr/include
LIBS += -L/usr/lib -lHepMC

#----------------------------------------------------------------------------------------------------

GCC = g++ -Wall -O3 -fPIC
ROOTCINT = rootcint

# modules for which dictionary shall be generated
dict_modules = Model IslamModel PPPModel BSWModel BHModel JenkovszkyModel ExpModel InterpolationModel CoulombInterference Constants

# modules with no dictionaries
modules = Math Generator

# executable files
exe_files = ElegentTest

#----------------------------------------------------------------------------------------------------

# directories to create
dirsToCreate = obj dict lib bin

dict_modules_obj = $(addprefix obj/, $(addsuffix .o, $(dict_modules)))
dict_modules_dict_src = $(addprefix dict/, $(addsuffix _dict.cc, $(dict_modules)))
dict_modules_dict_obj = $(addprefix dict/, $(addsuffix _dict.o, $(dict_modules)))

modules_obj = $(addprefix obj/, $(addsuffix .o, $(modules)))

exe_files_bin = $(addprefix bin/, $(exe_files))

#----------------------------------------------------------------------------------------------------

.PHONY: all dirs dicts libs exe info clean

all: libs exe

libs: dirs dicts lib/libElegent.so

exe: dirs libs $(exe_files_bin)

dirs: $(dirsToCreate)

$(dirsToCreate):
	mkdir -p $@

#dicts: $(dict_modules_dict_obj)

$(dict_modules_dict_obj) : dict/%_dict.o : dict/%_dict.cc
	$(GCC) $(INCS) -c $< -o $@

$(modules_dict_src) : dict/%_dict.cc : interface/%.h
	rm -f $@
	rm -f $(@:%.cc=%.h)
	$(ROOTCINT) $@ -c $<+

lib/libElegent.so: $(dict_modules_obj) $(modules_obj)
	@echo MAKING libElegent.so
	@$(GCC) -shared -olib/libElegent.so obj/*.o #dict/*.o

%.o : ../src/%.cc ../interface/%.h
	@echo COMPILING $<
	@$(GCC) $(INCS) -c $< -o $@

bin/% : src/%.cc lib/libElegent.so
	@echo BUILDING $@
	@$(GCC) $(INCS) $(LIBS) -Llib -lElegent $< -o $@

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
	rm -rf bin
