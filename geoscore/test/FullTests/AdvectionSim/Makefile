

.PHONY: dissolverTest ssppTest

all: dissolverTest ssppTest

dissolverTest:	
	(rm -f test_dissolver* ; \
	../../../src/GEOS.x -i DissolverTest.xml -m FlatFracture.geom; \
	../../../external_libs/bin/silodiff -q -l 1 test_dissolver_0100 benchmark/test_dissolver_0100 > dissolverTest.diff;);


ssppTest:	
	(rm -f test_sspp* ; \
	../../../src/GEOS.x -i SSParallelPlateTest.xml -m FlatFracture.geom; \
	../../../external_libs/bin/silodiff -q -l 1 test_sspp_0010 benchmark/test_sspp_0010 > ssppTest.diff;);

