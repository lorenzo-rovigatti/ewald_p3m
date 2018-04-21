#include <cstdlib>
#include <iostream>

#include "System.h"
#include "Ewald.h"

int main(int argc, char *argv[]) {
	srand48(1039098);

	System syst(10., 4);
	syst.print_conf("init.dat");
	Ewald ewald(syst);

	number energy = ewald.energy();
	std::cout << energy << " " << energy / syst.N() << std::endl;

	return 0;
}
