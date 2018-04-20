# MakeFile

# Set the Compiler
CC = g++

LIBDIR = ../../src/bdtlib/
XMLDIR = ../../src/tinyxml2/
CFLAGS = -I${LIBDIR} -I${XMLDIR} -c -O3 -g `root-config --cflags`
LIBFLAGS = `root-config --libs` -g -O3 -lXMLIO -lMLP -lMinuit 
DEPS = LoadEvents.hxx ${LIBDIR}SignificanceMetrics.hxx ${LIBDIR}Utilities.h
MAINPROGRAM = BasicTrainAndTest

CompileStuff: ${LIBDIR}CategoryReader.o ${LIBDIR}Node.o ${XMLDIR}tinyxml2.o ${LIBDIR}Tree.o ${MAINPROGRAM}.o
	$(CC) $(LIBFLAGS) ${LIBDIR}CategoryReader.o ${LIBDIR}Node.o ${XMLDIR}tinyxml2.o ${LIBDIR}Tree.o ${MAINPROGRAM}.o -o ${MAINPROGRAM}

${MAINPROGRAM}.o: ${MAINPROGRAM}.cxx ${DEPS}
	$(CC) $(CFLAGS) ${MAINPROGRAM}.cxx  

${LIBDIR}Node.o: ${LIBDIR}Node.cxx
	$(CC) $(CFLAGS) ${LIBDIR}Node.cxx -o ${LIBDIR}Node.o

${XMLDIR}tinyxml2.o: ${XMLDIR}tinyxml2.cpp
	$(CC) $(CFLAGS) ${XMLDIR}tinyxml2.cpp -o ${XMLDIR}tinyxml2.o

${LIBDIR}Tree.o: ${LIBDIR}Tree.cxx
	$(CC) $(CFLAGS) ${LIBDIR}Tree.cxx -o ${LIBDIR}Tree.o

${LIBDIR}CategoryReader.o: ${LIBDIR}CategoryReader.cxx
	$(CC) $(CFLAGS) ${LIBDIR}CategoryReader.cxx -o ${LIBDIR}CategoryReader.o

clean:
	rm *.o

cleanall:
	rm *.o
	rm ${LIBDIR}*.o
	rm ${XMLDIR}*.o
