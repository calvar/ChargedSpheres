BIN		:= md_ch_sphere

CUDA_PATH 	:= /usr/local/cuda-7.5
INCLUDES 	+= -I. -I$(CUDA_PATH)/include
LIBS		:= -L$(CUDA_PATH)/lib64

CXXFLAGS	:= -O3 -m64 #-g -fno-inline 
NVCCFLAGS	:= --default-stream per-thread #-g -G #allows per-thread default stream
LDFLAGS		:= -lgsl -lgslcblas -lgomp #-lcuda -lcudart

NVCC		:= nvcc
CXX		:= g++
LINKER		:= g++

C_SOURCES	:= main.cpp functions.cpp classes.cpp
CU_SOURCES	:= dev_functions.cu
HEADERS		:= functions.hpp classes.hpp dev_functions.hpp

C_OBJS		:= $(patsubst %.cpp, %.o, $(C_SOURCES))
CU_OBJS		:= $(patsubst %.cu, %.o, $(CU_SOURCES))


$(BIN): $(C_OBJS) $(HEADERS) #$(CU_OBJS)
	$(LINKER) -o $@ $(C_OBJS) $(LDFLAGS) $(INCLUDES) $(LIBS) #$(CU_OBJS)

$(C_OBJS): $(C_SOURCES) $(HEADERS)
	$(CXX) -c $(C_SOURCES) $(CXXFLAGS) $(INCLUDES) -fopenmp

#$(CU_OBJS): $(CU_SOURCES) $(HEADERS)
#	$(NVCC) -c $(CU_SOURCES) $(NVCCFLAGS) $(INCLUDES)

#sumProjection: sumProjection.cpp functions.hpp classes.hpp
#	$(CXX) $(OPTS) -w -o $@ sumProjection.cpp functions.cpp classes.cpp $(CXXFLAGS) -lgsl -lgslcblas -lgomp

graphic: graphic.cpp gl2ps.h
	$(CXX) $(OPTS) -w -o $@ graphic.cpp gl2ps.c -lGL -lGLU -lSDL

TimeCorrel: TimeCorrel.cpp 
	$(CXX) $(OPTS) -w -o $@ TimeCorrel.cpp $(CXXFLAGS)

integrate: integrate.cpp
	$(CXX) $(OPTS) -w -o $@ integrate.cpp $(CXXFLAGS)

trajectory: trajectory.cpp
	$(CXX) $(OPTS) -w -o $@ trajectory.cpp $(CXXFLAGS)

clean: 
	rm $(BIN) graphic *.o core*
