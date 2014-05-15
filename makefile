# MakeFile

# Set the Compiler
CC = g++

CFLAGS = -c `root-config --cflags`
LIBFLAGS = `root-config --libs` -lXMLIO -lMLP -lMinuit

FindBestParameters: Utilities.o Node.o Tree.o Forest.o FindBestParameters.o
	$(CC) $(LIBFLAGS) ./lib/Utilities.o ./lib/Node.o ./lib/Tree.o ./lib/Forest.o FindBestParameters.o -o FindBestParameters

FindBestParameters.o: FindBestParameters.cxx
	$(CC) $(CFLAGS) FindBestParameters.cxx  

Utilities.o: ./lib/Utilities.cxx
	$(CC) $(CFLAGS) ./lib/Utilities.cxx -o ./lib/Utilities.o

Node.o: ./lib/Node.cxx
	$(CC) $(CFLAGS) ./lib/Node.cxx -o ./lib/Node.o

Tree.o: ./lib/Tree.cxx
	$(CC) $(CFLAGS) ./lib/Tree.cxx -o ./lib/Tree.o

Forest.o: ./lib/Forest.cxx
	$(CC) $(CFLAGS) ./lib/Forest.cxx  -o ./lib/Forest.o


clean:
	rm *.o

cleanall:
	rm *.o
	rm ./lib/*.o
