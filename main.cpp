#include <cstdlib>
#include <iostream>

#include "System.h"
#include "Ewald.h"
#include "P3M.h"

int main(int argc, char *argv[]) {
	srand48(1039098);

	System syst(10., 10, 10.0, 10.0);
	syst.print_conf("init.dat");
	Ewald ewald(syst);
	P3M p3m(syst, 8, 1);

	std::cout << "EWALD" << std::endl;
	ewald.print_energy();
	std::cout << std::endl;

	std::cout << "P3M" << std::endl;
	p3m.print_energy();
	std::cout << std::endl;

	return 0;
}
