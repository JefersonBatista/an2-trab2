CC         = gcc
CFLAGS     = -g -c -lm -Wall
GCFLAGS    = -lm -Wall 
SOURCES    = ./CommonFiles/matrix.c        \
	     ./CommonFiles/linked_list.c   \
	     ./CommonFiles/graph.c         \
	     ./CommonFiles/ilup.c          \
	     ./CommonFiles/rcm.c           \
	     ./CommonFiles/gmres.c         \
	     program.c
OBJECTS    = $(SOURCES:.c=.o)
EXECUTABLE = program 

all: $(SOURCES) $(EXECUTABLE)
 
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(GCFLAGS) -g -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)



