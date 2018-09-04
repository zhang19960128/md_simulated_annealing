#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>
#include "interface.h"
int main(){
	std::fstream fs;
	fs.open("in.BTO",std::fstream::in);
	getparameter("pair_coeff",fs);
	return 0;
}
