
NODE_DATA_STRUCTURE:=0

###########################################################################
#
# Override default makefile settings
#
## Create or edit one of the following files:
##   Make.local.defs
##   Make.local.defs.(username) [e.g. Make.local.defs.benliu]
##   Make.local.defs.(host) [e.g. Make.local.defs.mcr]
##   Make.local.defs.(username).(host) [e.g. Make.local.defs.benliu.mcr]
## in the build directory to change the default compilation behavior
## `Make.local.defs.usertemplate' is provided as a template.

-include Make.local.defs
-include Make.local.defs.$(shell whoami)
###########################################################################

ifneq "$(wildcard src )" ""
$(info found src)
-include src/Makefile
CPP_FLAGS+=-Isrc -Isrc/Common
else
$(error didn't find src)
endif

ifneq "$(wildcard src_external )" ""
$(info found src_external)
  ifeq ($(MAKE_CORE_ONLY), 1)
    $(info ... but making core only)
  else
    -include src_external/Makefile
    CPP_FLAGS+=-Isrc_external -DSRC_EXTERNAL
  endif
else
  $(info didn't find src_external)
endif

ifneq "$(wildcard src_internal )" ""
$(info found src_internal)
  ifeq ($(MAKE_CORE_ONLY), 1)
    $(info ... but making core only)
  else
	ifeq ($(MAKE_EXTERNAL), 1)
$(info ... but making external only)
	else
		-include src_internal/Makefile
		CPP_FLAGS+=-Isrc_internal -DSRC_INTERNAL
	endif
  endif
else
$(info didn't find src_internal)
endif

OBJS = $(ARRAY_OBJS) $(BC_OBJS) $(CONTACT_OBJS) $(DATASTRUCTURES_OBJS) $(ELEMENT_OBJS) $(IC_OBJS) $(IO_OBJS) $(LINALG_OBJS) $(LOCAL_OBJS) $(MATERIAL_OBJS) $(MISC_OBJS) $(OBJMANAGER_OBJS) $(SOLVERS_OBJS) $(UTILITIES_OBJS) $(TOTALVIEW_OBJS)


ifeq ($(HAVE_MPI), 1)
	CPPC=$(CPP)
	OBJS += $(MPICOMM_OBJS)
else
	CPPC=$(SERIAL_CPP)
endif

CPP_FLAGS+=-I. -I$(EXTERNAL_LIBS)/include -DNODE_STRUCTURE=0 -DGPAC_MPI=$(HAVE_MPI) -DWRITE_FACENODE=$(WRITE_FACENODE)

LIBPATH_FLAGS+=-L$(EXTERNAL_LIBS)/lib
LD_FLAGS=-Wl,-rpath,$(EXTERNAL_LIBS)/lib

ifeq ($(HAVE_GMM), 1)
	CPP_FLAGS+=-DGPAC_GMM=$(HAVE_GMM) -DGD_MAT_LIB_CPP_INTERFACE
	CPP_FLAGS+=-I$(GMML_INC)/
        LNK_FLAGS+=-l$(GMML_LIB) 
	OBJS += $(GMM_OBJS)
	LIBPATH_FLAGS+=-L$(GMML_LIBPATH)
endif

ifeq ($(DEBUG), 1)
	CPP_FLAGS += -DRANGE_CHECKING=1  -DDEBUG=1
else
	CPP_FLAGS += -DRANGE_CHECKING=0 -DDEBUG=0
endif


ifeq ($(COMPILER), gcc)
	CPP_FLAGS +=  -Wall -fmessage-length=0
	ifeq ($(BUILD_TYPE), debug)
		CPP_FLAGS +=
	else
		CPP_FLAGS +=
	endif
	LIBPATH_FLAGS += -L/usr/lib64
	LIBPATH_FLAGS += -L/usr/lib/libblas	
endif
ifeq ($(COMPILER), icc)
	ifeq ($(BUILD_TYPE), debug)
		CPP_FLAGS +=
	else
		CPP_FLAGS +=
	endif
	LIBPATH_FLAGS += -L/usr/lib64
endif

# check for git
GIT_VERSION := $(shell git --version 2>/dev/null)

# get git data
ifdef GIT_VERSION
	REPO_VERSION := $(shell git describe --dirty --always)
	CPP_FLAGS += -DREPO_VERSION=\"$(REPO_VERSION)\"

	# may want to reduce length of repo history
	REPO_HISTORY += " $(shell echo \"\`git log --pretty=oneline\`\" ) "
	CPP_FLAGS += -DINCLUDE_VERSION_HISTORY=\"1\"


	ifneq "$(wildcard src_external )" ""
		ifneq ($(MAKE_CORE_ONLY), 1)
		  EXTERNAL_REPO_VERSION := $(shell (cd src_external/; git describe --dirty --always) )
		  CPP_FLAGS += -DEXTERNAL_REPO_VERSION=\"$(EXTERNAL_REPO_VERSION)\" 
		endif
        else
		EXTERNAL_REPO_VERSION := ""
	endif

	ifneq "$(wildcard src_internal )" ""
		ifneq ($(MAKE_CORE_ONLY), 1)
		  ifneq ($(MAKE_EXTERNAL), 1)
			INTERNAL_REPO_VERSION := $(shell (cd src_internal/; git describe --dirty --always) )
			CPP_FLAGS += -DINTERNAL_REPO_VERSION=\"$(INTERNAL_REPO_VERSION)\"
		  endif        
		endif        
        else
		INTERNAL_REPO_VERSION := ""
	endif
endif

all: GEOS.x repo_version external_repo_version internal_repo_version 

# the following forces recompilation of src common version if the git repository hashes change
.PHONY: force
repo_version: force 
	$(info Checking Core Repository Version)
	echo '$(REPO_VERSION)' | cmp -s - $@ || echo 'Updating core repository version'
	echo '$(REPO_VERSION)' | cmp -s - $@ || echo '$(REPO_VERSION)' > $@s
	$(info Recording Version History)
	echo "$(REPO_HISTORY)" > src/Common/VersionHistory.temp
	cat src/Common/VersionHistory.temp | sed "s/.*/ std::cout << \" & \" << std::endl; /" > src/Common/VersionHistory.temp.h

external_repo_version: force 
	$(info Checking External Repository Version)
	echo '$(EXTERNAL_REPO_VERSION)' | cmp -s - $@ || echo 'Updating external repository version'
	echo '$(EXTERNAL_REPO_VERSION)' | cmp -s - $@ || echo '$(EXTERNAL_REPO_VERSION)' > $@

internal_repo_version: force 
	$(info Checking Internal Repository Version)
	echo '$(INTERNAL_REPO_VERSION)' | cmp -s - $@ || echo 'Updating internal repository version'
	echo '$(INTERNAL_REPO_VERSION)' | cmp -s - $@ || echo '$(INTERNAL_REPO_VERSION)' > $@




GEOS.x: $(OBJS)
#	echo $(OBJS);
	$(CPPC) $(OBJS) $(LIBPATH_FLAGS)  $(LD_FLAGS) $(LNK_FLAGS) $(LNK_FLAGS_TRILINOS) -o GEOS.x && echo && echo SUCCESS;

Make.GEOS.config: $(OBJS) Makefile
	echo GEOS_OBJS=$(OBJS) > $@
	echo GEOS_LD_FLAGS=$(LNK_FLAGS) $(LNK_FLAGS_TRILINOS) $(LD_FLAGS) >> $@


# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)



# compile and generate dependency info;
# more complicated dependency computation, so all prereqs listed
# will also become command-less, prereq-less targets
#   sed:    strip the target (everything before colon)
#   sed:    remove any continuation backslashes
#   fmt -1: list words one per line
#   sed:    strip leading spaces
#   sed:    add trailing colons

%.o: %.cpp
ifeq ($(COMPILER), xlc)
	$(CPPC) -c $(CPP_FLAGS) $*.cpp -o $*.o
else
	$(CPPC) -c $(CPP_FLAGS) $*.cpp -o $*.o
	$(CPPC) -MM $(CPP_FLAGS) $*.cpp > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp
endif

#%.o: %.cpp
#	$(CPPC) -c $(CPP_FLAGS) $*.cpp -o $*.o
#	$(CPPC) -MM $(CPP_FLAGS) $*.cpp > $*.d
#	@mv -f $*.d $*.d.tmp
#	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
#	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
#	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
#	@rm -f $*.d.tmp

clean :
	rm -f GEOS.x $(OBJS)
	find . -type f -name "*.d" -print0 | xargs -0 rm -rf 


src/Common/Version.o: src/Common/Version.cpp  repo_version external_repo_version internal_repo_version
	$(CPPC) -c $(CPP_FLAGS) src/Common/Version.cpp -o src/Common/Version.o

