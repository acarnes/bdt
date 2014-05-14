# MakeFile

# Set the Compiler
CC = g++

CFLAGS = -c `root-config --cflags`
LIBFLAGS = `root-config --libs` -lXMLIO -lMLP -lMinuit

FindBestParameters: Utilities.o Node.o Tree.o Forest.o FindBestParameters.o
	$(CC) $(LIBFLAGS) Utilities.o Node.o Tree.o Forest.o FindBestParameters.o -o FindBestParameters

Utilities.o: Utilities.cxx
	$(CC) $(CFLAGS) Utilities.cxx

Node.o: Node.cxx
	$(CC) $(CFLAGS) Node.cxx

Tree.o: Tree.cxx
	$(CC) $(CFLAGS) Tree.cxx  

Forest.o: Forest.cxx
	$(CC) $(CFLAGS) Forest.cxx  

FindBestParameters.o: FindBestParameters.cxx
	$(CC) $(CFLAGS) FindBestParameters.cxx  

clean:
	rm *.o
