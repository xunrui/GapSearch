
TARGETS = gapsearch.exe \

OBJECTS = 	

CC 	= gcc
CXX = g++

OPT =  -O2
# -O3

DEBUG =  
# -g 
# -static 
WARNINGS = -Wall


CFLAGS 		= $(OPT) $(DEBUG) $(WARNINGS)
CXXFLAGS 	= $(OPT) $(DEBUG) $(WARNINGS)
LFLAGS 		= $(OPT) $(DEBUG) $(WARNINGS)

LIBS 	 = -lglpk

INCLUDES =  -I../nauty24r2 


TREESEARCHOBJS	= 
				  
NAUTYOBJS     	= ../nauty24r2/nauty.o 			\
				  ../nauty24r2/nausparse.o		\
				  ../nauty24r2/gtools.o			\
				  ../nauty24r2/nautil.o			\
				  ../nauty24r2/naugraph.o
					
LIBOBJS			= $(TREESEARCHOBJS) $(NAUTYOBJS)
			
			
.SUFFIXES: .c .cpp .o .obj .exe 

all: $(OBJECTS) $(TESTS) $(TARGETS)


# The default object compiler
.c.o: $<
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -c $< -o $@
        
.cpp.o: $<
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS) -c $< -o $@
        
.cpp.exe: $< $(OBJECTS)  
	$(CXX) $(LFLAGS)			\
        	$(INCLUDES)	$(DEBUG)			\
        	$(LIBOBJS) $(LIBS)				\
        	`cat $@.objs`           		\
            $< -o $@
        
.c.exe: $< $(COBJECTS)
	$(CC) 	$(LFLAGS)			    \
        	$(INCLUDES)				\
        	$(NAUTYOBJS)  $(COBJECTS) $(LIBS)		\
            $< -o $@
        
clean:
	rm $(OBJECTS) $(TARGETS) $(TESTS)
	
cleanexe:
	rm $(TARGETS)

clexe:
	rm $(TARGETS)
