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

  // Write to screen how the experiment is configured
  std::cout << "==============================" << std::endl;
  std::cout << "|    How am I configured?    |" << std::endl;
  std::cout << "==============================" << std::endl;
  config.Write(std::cout);
  std::cout << "==============================\n" << std::endl;

  emp::Random rnd(config.SEED());

  if (config.PROBLEM() == (uint32_t)PROBLEM_TYPE::REAL_VALUE) {
    using ORG_TYPE = rv_t;
    EcologyWorld<ORG_TYPE> world(rnd);
    world.Setup(config);   
    if (config.MODE() == 0) {
      world.Run(); 
    } else {
      world.TimeCompetition();
    }
  } else if (config.PROBLEM() == (uint32_t)PROBLEM_TYPE::NK) {
    using ORG_TYPE = bit_t;
    EcologyWorld<ORG_TYPE> world(rnd);
    world.Setup(config);
    if (config.MODE() == 0) {
      world.Run(); 
    } else {
      world.TimeCompetition();
    }
  } else if ((config.PROBLEM() == (uint32_t)PROBLEM_TYPE::LOGIC_9) || 
             (config.PROBLEM() == (uint32_t)PROBLEM_TYPE::PROGRAM_SYNTHESIS)) {
    using ORG_TYPE = gp_t;
    EcologyWorld<ORG_TYPE> world(rnd);
    world.Setup(config);
    if (config.MODE() == 0) {
      world.Run(); 
    } else {
      world.TimeCompetition();
    }
  } else if (config.PROBLEM() == (uint32_t)PROBLEM_TYPE::SORTING_NETWORK) {
    using ORG_TYPE = sorting_t;
    EcologyWorld<ORG_TYPE> world(rnd);
    world.Setup(config);
    if (config.MODE() == 0) {
      world.Run(); 
    } else {
      world.TimeCompetition();
    }
  }
}
