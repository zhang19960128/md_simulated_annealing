#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>
#include "readion.h"
#include <ctime>
#include "interface.h"
#include "readpara.h"
#include "penalty.h"
int main(){
	 readPT("control.PT");
	 readvmmap("param.map");
	 readbound("param.vmbound");
	 int size_box;
	 readion(control::ionfile[0],40,size_box,8);
	 return 0;
	 
}
