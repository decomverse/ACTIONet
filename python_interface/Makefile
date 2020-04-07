LIBNAME:=libactionet.a

CXX=dpcpp
CXXFLAGS=-O2 -std=c++17 -fsycl -pthread -fPIC -w -DMKL_ILP64 
EXTFLAGS=-fsycl-unnamed-lambda


#LINALG=-lopenblas -llapack
#LINALG=-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl	
LINALG=-L${MKLROOT}/lib/intel64 -lmkl_sycl -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lsycl -lOpenCL -lpthread -lm -ldl
LIB_FLAGS+=-lstdc++ ${LINALG}
INCLUDE=-I./include/ -I./include/arma/ -I./include/ACTIONet/ -I./include/ACTIONet/SPAMS -I./include/ACTIONet/hnsw -I${MKLROOT}/include

SRC=$(shell find src/ACTIONet -type f -name "*.cc")
OBJ=$(SRC:.cc=.o)

PROGRAM=run_test
	
all: $(PROGRAM) message
	
src/%.o: src/%.cc
	$(CXX) $(CXXFLAGS) ${INCLUDE} -c -o $@  $<


run_test: $(OBJ) src/run_test.o
	$(CXX) $(CXXFLAGS)  -o $@ src/run_test.o -L. $(OBJ) $(LIB_FLAGS)	

.PHONY: clean ar

clean:
	rm -f $(PROGRAM) $(OBJ) $(LIBNAME) 

ar:
	make clean
	tar -czvf ../$(PROGRAM)"(`date`)".tar.gz *

message:
	@echo "Executables: $(PROGRAM) have been created"	
