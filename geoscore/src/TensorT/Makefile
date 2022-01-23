CPPCP = g++
LNK = g++


#CPPFLAGS = -O3 -Wall -funroll-loops  -mfpmath=sse -msse3 -march=core2
CPPFLAGS = -O3 -Wall -mfpmath=sse -msse3 -march=core2
#CPPFLAGS = -O0 -g3 
LNKFLAGS =

OBJS = compund_operation_benchmark.o \
nodal_update_benchmark.o

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

all: nodal_update_benchmark.x compund_operation_benchmark.x

nodal_update_benchmark.x: nodal_update_benchmark.o
	$(LNK) $(LNKFLAGS) nodal_update_benchmark.o -o nodal_update_benchmark.x

compund_operation_benchmark.x: compund_operation_benchmark.o
	$(LNK) $(LNKFLAGS) compund_operation_benchmark.o -o compund_operation_benchmark.x




# compile and generate dependency info;
# more complicated dependency computation, so all prereqs listed
# will also become command-less, prereq-less targets
#   sed:    strip the target (everything before colon)
#   sed:    remove any continuation backslashes
#   fmt -1: list words one per line
#   sed:    strip leading spaces
#   sed:    add trailing colons

%.o: %.cpp
	$(CPPCP) -c $(CPPFLAGS) $*.cpp -o $*.o
	$(CPPCP) -MM $(CPPFLAGS) $*.cpp > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp


# pull in dependency info for *existing* .o files
#-include $(OBJS:.o=.P)
#%.o : %.cpp
#	$(CPPCP) $(CPPFLAGS) -MD -o $@ $<
#	@cp $*.d $*.P; \
#	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
#	rm -f $*.d


clean :
	rm *.x *.o *.d
