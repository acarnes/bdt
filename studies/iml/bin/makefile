# MakeFile

# Set the Compiler
CC = g++

LIBDIR = ../../../src/bdtlib/
XMLDIR = ../../../src/tinyxml2/
CFLAGS = -I${LIBDIR} -I${XMLDIR} -c -O3 -g `root-config --cflags`
LIBFLAGS = `root-config --libs` -g -O3 -lXMLIO -lMLP -lMinuit 
DEPS = LoadSaveEvents.h Functions.h ${XMLDIR}tinyxml2.h
MAINPROGRAM = BasicTrainAndTest

CompileStuff: ${LIBDIR}Utilities.o ${LIBDIR}Node.o ${XMLDIR}tinyxml2.o ${LIBDIR}Tree.o ${LIBDIR}Forest.o ${MAINPROGRAM}.o
	$(CC) $(LIBFLAGS) ${LIBDIR}Utilities.o ${LIBDIR}Node.o ${XMLDIR}tinyxml2.o ${LIBDIR}Tree.o ${LIBDIR}Forest.o ${MAINPROGRAM}.o -o ${MAINPROGRAM}

${MAINPROGRAM}.o: ${MAINPROGRAM}.cxx ${DEPS}
	$(CC) $(CFLAGS) ${MAINPROGRAM}.cxx  

${LIBDIR}Utilities.o: ${LIBDIR}Utilities.cxx
	$(CC) $(CFLAGS) ${LIBDIR}Utilities.cxx -o ${LIBDIR}Utilities.o

${LIBDIR}Node.o: ${LIBDIR}Node.cxx
	$(CC) $(CFLAGS) ${LIBDIR}Node.cxx -o ${LIBDIR}Node.o

${XMLDIR}tinyxml2.o: ${XMLDIR}tinyxml2.cpp
	$(CC) $(CFLAGS) ${XMLDIR}tinyxml2.cpp -o ${XMLDIR}tinyxml2.o

${LIBDIR}Tree.o: ${LIBDIR}Tree.cxx ${XMLDIR}libtinyxml2.a
	$(CC) $(CFLAGS) ${LIBDIR}Tree.cxx -o ${LIBDIR}Tree.o

${LIBDIR}Forest.o: ${LIBDIR}Forest.cxx ${LIBDIR}LossFunctions.h
	$(CC) $(CFLAGS) ${LIBDIR}Forest.cxx  -o ${LIBDIR}Forest.o

clean:
	rm *.o

cleanall:
	rm *.o
	rm ${LIBDIR}*.o
	rm ${XMLDIR}*.o
