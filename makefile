# MakeFile

# Set the Compiler
CC = g++

CFLAGS = -c `root-config --cflags`
LIBFLAGS = `root-config --libs` -lXMLIO -lMLP -lMinuit

Forest: Node.o Tree.o DrawPlots.o Forest.o
	$(CC) $(LIBFLAGS) Node.o Tree.o DrawPlots.o Forest.o -o Forest

Node.o: Node.cxx
	$(CC) $(CFLAGS) Node.cxx

Tree.o: Tree.cxx
	$(CC) $(CFLAGS) Tree.cxx  

DrawPlots.o: DrawPlots.cxx
	$(CC) $(CFLAGS) DrawPlots.cxx  

Forest.o: Forest.cxx
	$(CC) $(CFLAGS) Forest.cxx  

clean:
	rm *.o
