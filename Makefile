MPI_CC = mpicxx
LD = mpicxx
CCFLAGS = -std=c++11 -c -Wall 
LDFLAGS = -lfftw_mpi -lfftw

TARGET = fil.e
OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))

SHELL := /bin/bash



all: $(TARGET)

$(TARGET): $(OBJECTS)
	. /usr/share/Modules/init/bash; \
	module load openmpi/1.8.4-icc; \
	module load intel/15.0.090; \
	$(LD) -o $@ $^ $(LDFLAGS) -I/home/users/dergachev/local/include -L/home/users/dergachev/local/lib

%.o: %.cpp
	. /usr/share/Modules/init/bash; \
	module load openmpi/1.8.4-icc; \
	module load intel/15.0.090; \
	$(MPI_CC) $(CCFLAGS) $^ -o $@ -I/home/users/dergachev/local/include
	
clean:
	rm $(TARGET) $(OBJECTS)