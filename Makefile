CXX=icc
CXXFLAGS= -std=c++11
DEPS= interface.h image.h atom.h
LIBPATH =
CXXFLAGS +=
CXXFLAGS +=-I./include
VPATH=./src
sa.x: interface.o atom.o main.o ewald_jiahaoz.o lj12.o bv.o bvv.o image.o stress.o readion.o
	mkdir -p obj bin
	$(CXX) -o sa.x $(LIBPATH) interface.o atom.o main.o ewald_jiahaoz.o lj12.o bv.o bvv.o image.o stress.o readion.o
	mv *.o obj
	mv sa.x bin
%.o: %.c $(DEPS)
	$(CXX)  $(CXXFLAGS) -c -o $@ $^
clean:
	rm -rf *.o obj bin
