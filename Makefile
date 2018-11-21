CXX=g++
CXXFLAGS= -std=c++11
DEPS= interface.h image.h atom.h readion.h readpara.h sa.h penalty.h simann.h
LIBPATH =
CXXFLAGS +=
CXXFLAGS +=-I./include
vpath %.cpp src
vpath %.h include
vpath %.o obj
sa.x: interface.o atom.o ewald_jiahaoz.o lj12.o bv.o bvv.o image.o stress.o readion.o readpara.o main.o sa.o simann.o penalty.o
	mkdir -p obj bin
	$(CXX) -o sa.x $(LIBPATH) interface.o atom.o main.o ewald_jiahaoz.o lj12.o bv.o bvv.o image.o stress.o readion.o readpara.o sa.o simann.o penalty.o
	mv *.o bin/
%.o: %.c $(DEPS)
	$(CXX)  $(CXXFLAGS) -c -o $@ $^
clean:
	rm -rf *.o obj bin
