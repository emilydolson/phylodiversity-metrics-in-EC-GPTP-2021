// This is the main function for the NATIVE version of this project.

#include <iostream>

#include "../ecology_world.h"

#include "base/vector.h"
#include "config/command_line.h"

int main(int argc, char* argv[])
{
  emp::vector<std::string> args = emp::cl::args_to_strings(argc, argv);

  EcologyWorld<int> world;

  std::cout << "Hello World!" << std::endl;
}
