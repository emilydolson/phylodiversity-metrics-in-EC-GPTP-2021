// This is the main function for the NATIVE version of this project.

#include <iostream>

#include "../ecology_world.h"

#include "base/vector.h"
#include "config/command_line.h"

int main(int argc, char* argv[])
{
  EcologyConfig config;
  auto args = emp::cl::ArgManager(argc, argv);
  if (args.ProcessConfigOptions(config, std::cout, "EcologyConfig.cfg", "Ecology-macros.h") == false) exit(0);
  if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

  emp::Random rnd(config.SEED());

  if (config.PROBLEM() == (uint32_t)PROBLEM_TYPE::REAL_VALUE) {
    using ORG_TYPE = emp::vector<double>;
    EcologyWorld<ORG_TYPE> world(rnd);
    world.Setup(config);
    world.Run();
  } else if (config.PROBLEM() == (uint32_t)PROBLEM_TYPE::NK) {
    using ORG_TYPE = emp::BitVector;
    EcologyWorld<ORG_TYPE> world(rnd);
    world.Setup(config);
    world.Run();
  }

}
