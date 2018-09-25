CXX=icc
CXXFLAGS= -std=c++11
DEPS= interface.h image.h atom.h
LIBPATH =-L/opt/gcc/7.3.0/snos/lib64
CXXFLAGS += -I/opt/gcc/7.3.0/snos/include
%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)
sa.x: interface.o atom.o main.o ewald_jiahaoz.o lj12.o bv.o bvv.o image.o
	$(CXX) -o sa.x $(LIBPATH) interface.o atom.o main.o ewald_jiahaoz.o lj12.o bv.o bvv.o image.o
clean:
	rm -rf *.o
