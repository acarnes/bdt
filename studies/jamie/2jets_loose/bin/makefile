# MakeFile

# Set the Compiler
CC = g++

LIBFLAGS = `root-config --libs` -lXMLIO -lMLP -lMinuit
LIBDIR = ../../../../lib/
CFLAGS = -I${LIBDIR} -c `root-config --cflags`
DEPS = LoadSaveEvents.h Functions.h
MAINPROGRAM = FindBestParameters

CompileStuff: Utilities.o Node.o Tree.o Forest.o ${MAINPROGRAM}.o
	$(CC) $(LIBFLAGS) ${LIBDIR}Utilities.o ${LIBDIR}Node.o ${LIBDIR}Tree.o ${LIBDIR}Forest.o ${MAINPROGRAM}.o -o ${MAINPROGRAM}

${MAINPROGRAM}.o: ${MAINPROGRAM}.cxx ${DEPS}
	$(CC) $(CFLAGS) ${MAINPROGRAM}.cxx  

Utilities.o: ${LIBDIR}Utilities.cxx
	$(CC) $(CFLAGS) ${LIBDIR}Utilities.cxx -o ${LIBDIR}Utilities.o

Node.o: ${LIBDIR}Node.cxx
	$(CC) $(CFLAGS) ${LIBDIR}Node.cxx -o ${LIBDIR}Node.o

Tree.o: ${LIBDIR}Tree.cxx
	$(CC) $(CFLAGS) ${LIBDIR}Tree.cxx -o ${LIBDIR}Tree.o

Forest.o: ${LIBDIR}Forest.cxx
	$(CC) $(CFLAGS) ${LIBDIR}Forest.cxx  -o ${LIBDIR}Forest.o


clean:
	rm *.o

cleanall:
	rm *.o
	rm ${LIBDIR}*.o
