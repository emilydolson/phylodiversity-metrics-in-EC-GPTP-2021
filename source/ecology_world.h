//  Copyright (C) Emily Dolson 2018.
//  Released under the MIT Software license; see doc/LICENSE
//
//
//  This file contains the bulk of the code for studying open-ended evolution in NK Landscapes

#include <iostream>

#include "config/ArgManager.h"
#include "Evolve/NK.h"
#include "Evolve/World.h"
#include "Evolve/Resource.h"
#include "tools/BitVector.h"
#include "tools/Random.h"
#include "tools/sequence_utils.h"
#include "Evolve/OEE.h"
#include "tools/vector_utils.h"
#include "tools/Graph.h"

#include "cec2013.h"

EMP_BUILD_CONFIG( EcologyConfig,
  GROUP(MAIN, "Global settings"),
  VALUE(SEED, int, 0, "Random number seed (0 for based on time)"),
  VALUE(POP_SIZE, uint32_t, 1000, "Number of organisms in the popoulation."),
  VALUE(MAX_GENS, uint32_t, 2000, "How many generations should we process?"),
  VALUE(MUT_RATE, double, .005, "Probability of each site being mutated. For real-valued problems, the standard deviation of of the distribution from which mutations are pulled."),
  VALUE(PROBLEM, uint32_t, 0, "Which problem to use? 0 = NK, 1 = Program synthesis, 2 = Real-valued, 3 = Sorting network, 4 = Logic-9"),

  GROUP(NK, "Settings for NK landscape"),
  VALUE(K, uint32_t, 10, "Level of epistasis in the NK model"),
  VALUE(N, uint32_t, 200, "Number of bits in each organisms (must be > K)"),

  GROUP(REAL_VALUED, "Settings for real-valued optimzation problems"),
  VALUE(FUNCTION_NUMBER, uint32_t, 0, "Problem to use"),
  VALUE(DIMS, uint32_t, 200, "Number of dimensions in orgs"),

  GROUP(TESTCASES, "Settings for problems that use testcases"),
  VALUE(TRAINSET_FILE, std::string, "testcases/count-odds.csv", "Which set of testcases should we use for training?"),  
  VALUE(TESTSET_FILE, std::string, "testcases/count-odds.csv", "Which set of testcases should we use for evaluation?"),
  VALUE(N_TEST_CASES, uint32_t, 200, "How many test cases to use"),  

  GROUP(SELECTION_METHODS, "Settings related to selection"),
  VALUE(SELECTION, uint32_t, 0, "Selection method. 0 = Tournament, 1 = fitness sharing, 2 = lexicase, 3 = Eco-EA, 4 = Random"),
  VALUE(TOURNAMENT_SIZE, int, 2, "For tournament selection, number of individuals to include in tournament"),
  VALUE(SHARING_ALPHA, double, 1, "Alpha for fitness sharing (controls shape of sharing function)"),
  VALUE(SHARING_THRESHOLD, double, 2, "For fitness sharing - how similar do two individuals have to be to compete?"),
  VALUE(RESOURCE_SELECT_RES_AMOUNT, double, 50.0, "Initial resource amount (for all resources)"),
  VALUE(RESOURCE_SELECT_RES_INFLOW, double, 50.0, "Resource in-flow (amount)"),
  VALUE(RESOURCE_SELECT_RES_OUTFLOW, double, 0.05, "Resource out-flow (rate)"),
  VALUE(RESOURCE_SELECT_FRAC, double, 0.0025, "Fraction of resource consumed."),
  VALUE(RESOURCE_SELECT_MAX_BONUS, double, 5.0, "What's the max bonus someone can get for consuming a resource?"),
  VALUE(RESOURCE_SELECT_COST, double, 0.0, "Cost of using a resource?"),

  GROUP(OPEN_ENDED_EVOLUTION, "Settings related to tracking MODES metrics"),
  VALUE(MODES_RESOLUTION, int, 1, "How often should MODES metrics be calculated?"),
  VALUE(FILTER_LENGTH, int, 1000, "How many generations should we use for the persistence filter?"),

  GROUP(DATA_RESOLUTION, "How often should data get printed?"),
  VALUE(FAST_DATA_RES, int, 10, "How often should easy to calculate metrics (e.g. fitness) be calculated?"),
  VALUE(ECOLOGY_DATA_RES, int, 100, "How often should ecological interactions (expensive) be calculated?"),
  VALUE(ECOLOGY_DATA_N_SIMULATIONS, int, 500, "How many selection events be simulated to calculate interactions?")
);

// Special version for bitstrings
// The reason the org can't be const is that it needs to get plugged into the fitness function, 
// which may not be const
emp::vector<int> Skeletonize(emp::BitVector & org, std::function<double(emp::BitVector&)> fun_calc_fitness) {
    emp_assert(org.size() > 0, "Empty org passed to skeletonize");

    emp::vector<int> skeleton;
    double fitness = fun_calc_fitness(org);
    emp::BitVector test_org = emp::BitVector(org);

    for (size_t i = 0; i < org.size(); i++) {
        test_org[i] = !test_org[i]; // For bit strings we just flip the bit
        double new_fitness = fun_calc_fitness(test_org);
        if (new_fitness < fitness) {
            skeleton.push_back(org[i]);
        } else {
            skeleton.push_back(-1);
        }
        test_org[i] = (int)org[i];
    }

    return skeleton;
}

template <typename ORG>
class EcologyWorld : public emp::World<ORG> {

    using typename emp::World<ORG>::fun_calc_fitness_t;
    using emp::World<ORG>::GetSize;
    using emp::World<ORG>::update;
    using emp::World<ORG>::random_ptr;
    using emp::World<ORG>::SetupFitnessFile;
    using emp::World<ORG>::SetupSystematicsFile;
    using emp::World<ORG>::SetupPopulationFile;
    using emp::World<ORG>::SetPopStruct_Mixed;
    using emp::World<ORG>::SetSynchronousSystematics;
    using emp::World<ORG>::pop;
    using emp::World<ORG>::fun_calc_fitness;
    using emp::World<ORG>::AddDataNode;
    using emp::World<ORG>::GetDataNode;
    using emp::World<ORG>::AddDataFile;
    using emp::World<ORG>::SetAutoMutate;
    using emp::World<ORG>::SetCache;
    using emp::World<ORG>::IsOccupied;
    using emp::World<ORG>::CalcFitnessID;
    using emp::World<ORG>::GetOrg;
    using emp::World<ORG>::GetRandomOrgID;
    using emp::World<ORG>::GetRandom;
    using emp::World<ORG>::Update;
    using emp::World<ORG>::GetGenomeAt;

    enum class SELECTION_METHOD { TOURNAMENT=0, SHARING=1, LEXICASE=2, ECOEA=3, RANDOM=4 };
    enum class PROBLEM_TYPE { NK=0, PROGRAM_SYNTHESIS=1, REAL_VALUE=2, SORTING_NETWORK=3, LOGIC_9=4 };

    emp::Signal<void(void)> evaluate_performance_sig;    ///< Triggered when we want to evaluate performance

    uint32_t N;
    uint32_t K;
    uint32_t POP_SIZE;
    uint32_t MAX_GENS;
    uint32_t SELECTION;
    uint32_t CHANGE_TYPE;
    uint32_t CHANGE_RATE;
    uint32_t PROBLEM;
    uint32_t N_TEST_CASES;
    double MUT_RATE;
    int TOURNAMENT_SIZE;
    int MODES_RESOLUTION;
    int FILTER_LENGTH;
    double SHARING_ALPHA;
    double SHARING_THRESHOLD;
    int FAST_DATA_RES;
    int ECOLOGY_DATA_RES;
    int ECOLOGY_DATA_N_SIMULATIONS;

    double RESOURCE_SELECT_RES_AMOUNT;
    double RESOURCE_SELECT_RES_INFLOW;
    double RESOURCE_SELECT_RES_OUTFLOW;
    double RESOURCE_SELECT_FRAC;
    double RESOURCE_SELECT_MAX_BONUS;
    double RESOURCE_SELECT_COST;

    emp::NKLandscape landscape;
 
    // Function to determine actual performance on problem
    fun_calc_fitness_t evaluation_fun;

    // Set of fitness functions for lexicase and Eco-EA
    emp::vector<fun_calc_fitness_t> fit_set;
    emp::vector<emp::Resource> resources;

    emp::Ptr<emp::OEETracker<emp::Systematics<ORG, ORG>, emp::vector<int>>> oee;


    // Data tracking
    emp::DataNode<double, emp::data::Histogram, emp::data::Stats> interaction_strengths;
    emp::DataNode<int, emp::data::Histogram, emp::data::Stats> node_in_degrees;
    emp::DataNode<int, emp::data::Histogram, emp::data::Stats> node_out_degrees;

    public:

    EcologyWorld() {;}
    EcologyWorld(emp::Random & rnd) : emp::World<ORG>(rnd) {;}
    ~EcologyWorld(){if (oee){oee.Delete();}}

    void Setup(EcologyConfig & config) {
        N = config.N();
        K = config.K();
        POP_SIZE = config.POP_SIZE();
        MAX_GENS = config.MAX_GENS();
        MUT_RATE = config.MUT_RATE();
        SELECTION = config.SELECTION();
        TOURNAMENT_SIZE = config.TOURNAMENT_SIZE();
        MODES_RESOLUTION = config.MODES_RESOLUTION();
        FAST_DATA_RES = config.FAST_DATA_RES();
        ECOLOGY_DATA_RES = config.ECOLOGY_DATA_RES();
        ECOLOGY_DATA_N_SIMULATIONS = config.ECOLOGY_DATA_N_SIMULATIONS();
        FILTER_LENGTH = config.FILTER_LENGTH();
        SHARING_THRESHOLD = config.SHARING_THRESHOLD();
        SHARING_ALPHA = config.SHARING_ALPHA();
        PROBLEM = config.PROBLEM();
        N_TEST_CASES = config.N_TEST_CASES();

        RESOURCE_SELECT_RES_AMOUNT = config.RESOURCE_SELECT_RES_AMOUNT();
        RESOURCE_SELECT_RES_INFLOW = config.RESOURCE_SELECT_RES_INFLOW();
        RESOURCE_SELECT_RES_OUTFLOW = config.RESOURCE_SELECT_RES_OUTFLOW();
        RESOURCE_SELECT_FRAC = config.RESOURCE_SELECT_FRAC();
        RESOURCE_SELECT_MAX_BONUS = config.RESOURCE_SELECT_MAX_BONUS();
        RESOURCE_SELECT_COST = config.RESOURCE_SELECT_COST();

        switch (PROBLEM){
            case (uint32_t) PROBLEM_TYPE::NK:
                SetupNK();
                break;

            case (uint32_t) PROBLEM_TYPE::PROGRAM_SYNTHESIS:
                // TODO
                break;            

            case (uint32_t) PROBLEM_TYPE::REAL_VALUE:
                // TODO
                break;            

            case (uint32_t) PROBLEM_TYPE::SORTING_NETWORKS:
                // TODO
                break;            

            case (uint32_t) PROBLEM_TYPE::LOGIC_9:
                // TODO
                break;            

            default:
                emp_assert(false && "Invalid problem type", PROBLEM);
                exit(1);
                break;
        }
 
        switch (SELECTION) {
            case ( (uint32_t) SELECTION_METHOD::ECOEA):
                resources.resize(0);
                for (size_t i=0; i<N_TEST_CASES; i++) {
                    resources.push_back(emp::Resource(RESOURCE_SELECT_RES_INFLOW, RESOURCE_SELECT_RES_INFLOW, .01));
                }
                break;

            default:
                break;

        }

        emp::Ptr<emp::Systematics<ORG, ORG> > sys;
        sys.New([](const ORG & o){return o;});
        oee.New(sys, [this](ORG & org){return Skeletonize(org, fun_calc_fitness);}, [](const emp::vector<int> & org){
            return org.size() - emp::Count(org, -1);
        });
        oee->SetResolution(MODES_RESOLUTION);
        oee->SetGenerationInterval(FILTER_LENGTH);
        AddSystematics(sys);
        OnUpdate([this](int ud){oee->Update(ud);});

        SetupFitnessFile().SetTimingRepeat(10);
        SetupSystematicsFile().SetTimingRepeat(10);
        SetupPopulationFile().SetTimingRepeat(10);
        SetPopStruct_Mixed(true);
        SetSynchronousSystematics(true);

        AddDataNode("performance");
        evaluate_performance_sig.AddAction([this](){
            auto perf_node = GetDataNode("performance");
            perf_node.Reset();
            for (ORG & org : pop) {
                perf_node.Add(evaluation_fun(org));
            }
        });

        emp::DataFile performance_file;
        performance_file.AddVar(update, "update", "Update");
        performance_file.AddStats(GetDataNode("performance"), "performance", "How well are programs actually solving the problem?");
        performance_file.PrintHeaderKeys();
        performance_file.SetTimingRepeat(FAST_DATA_RES);
        AddDataFile(performance_file);


        emp::DataFile oee_file;
        oee_file.AddVar(update, "generation", "Generation");
        oee_file.AddCurrent(*oee->GetDataNode("change"), "change", "change potential");
        oee_file.AddCurrent(*oee->GetDataNode("novelty"), "novelty", "novelty potential");
        oee_file.AddCurrent(*oee->GetDataNode("diversity"), "ecology", "ecology potential");
        oee_file.AddCurrent(*oee->GetDataNode("complexity"), "complexity", "complexity potential");
        oee_file.PrintHeaderKeys();
        oee_file.SetTimingRepeat(MODES_RESOLUTION);
        AddDataFile(oee_file);

        interaction_strengths.SetupBins(-1.0, 1.0, 10);
        node_out_degrees.SetupBins(0, POP_SIZE, 20);
        node_in_degrees.SetupBins(0, POP_SIZE, 20);

        emp::DataFile ecology_file;
        ecology_file.AddVar(update, "generation", "Generation");
        ecology_file.AddStats(interaction_strengths, "interaction_strength", "interaction strength");
        ecology_file.AddAllHistBins(interaction_strengths, "interaction_strength", "interaction strength");
        ecology_file.AddStats(node_out_degrees, "interaction_node_out_degree", "interaction node out degree"); 
        ecology_file.AddAllHistBins(node_out_degrees, "interaction_node_out_degree", "interaction node out degree"); 
        ecology_file.AddStats(node_in_degrees, "interaction_node_in_degree", "interaction node in degree");
        ecology_file.AddAllHistBins(node_in_degrees, "interaction_node_in_degree", "interaction node in degree");
        ecology_file.PrintHeaderKeys();
        ecology_file.SetTimingRepeat(MODES_RESOLUTION);
        AddDataFile(ecology_file);

        SetAutoMutate();
        if (SELECTION == (uint32_t)SELECTION_METHOD::TOURNAMENT) {
            SetCache();
        }

    }

    void CalcInteractionNetwork() {
        // Get estimate of current fitnesses


    }


    void SetupNK() {
        // Create landscape based on correct random seed
        landscape = emp::NKLandscape(N, K, *random_ptr);

        // Set-up fitness function (queries NK Landscape)
        evaluation_fun =
            [this](ORG & org){ return landscape.GetFitness(org); };

        if (SELECTION == (uint32_t) SELECTION_METHOD::SHARING) {
            SetSharedFitFun(fun_calc_fitness, [](ORG & org1, ORG & org2){return calc_hamming_distance(org1, org2);}, SHARING_THRESHOLD, SHARING_ALPHA);
        } else {
            SetFitFun(fun_calc_fitness);
        }

        // Make a fitness function for each gene (set of sites that determine a fitness contribution)
        // in the bitstring
        if (SELECTION == (uint32_t) SELECTION_METHOD::LEXICASE || SELECTION == (uint32_t) SELECTION_METHOD::ECOEA ) {
            for (size_t gene = 0; gene < N; gene++) {
                fit_set.push_back([this, gene](ORG & org){return landscape.GetFitness(gene, org);});
            }
        }

        // Build a random initial population
        for (uint32_t i = 0; i < POP_SIZE; i++) {
            ORG next_org(N);
            for (uint32_t j = 0; j < N; j++) next_org[j] = random_ptr->P(0.5);
            Inject(next_org);
        }

        // Setup the mutation function.
        std::function<size_t(ORG &, emp::Random &)> mut_fun =
            [this](ORG & org, emp::Random & random) {
                size_t num_muts = 0;
                for (uint32_t m = 0; m < N; m++) {
                    if (random_ptr->P(MUT_RATE)) {
                        org[m] = random_ptr->P(.5); // Randomly assign 0 or 1 
                        num_muts++;
                    }
                }
                return num_muts;
            };
        SetMutFun( mut_fun );

    }

    void RunStep() {

        std::cout << update << std::endl;
        switch(SELECTION) {
            case (uint32_t)SELECTION_METHOD::TOURNAMENT :
               emp::TournamentSelect(*this, TOURNAMENT_SIZE, POP_SIZE);
               break;

            case (uint32_t)SELECTION_METHOD::SHARING : // Sharing is handled in the setting of the fitness function
               emp::TournamentSelect(*this, TOURNAMENT_SIZE, POP_SIZE);
               break;

            case (uint32_t)SELECTION_METHOD::ECOEA :
                emp::ResourceSelect(*this, fit_set, resources, TOURNAMENT_SIZE, POP_SIZE, RESOURCE_SELECT_FRAC, RESOURCE_SELECT_MAX_BONUS,RESOURCE_SELECT_COST, true);
                break;

            case (uint32_t)SELECTION_METHOD::RANDOM :
                emp::RandomSelect(*this, POP_SIZE);
                break;

            case (uint32_t)SELECTION_METHOD::LEXICASE :
                emp::LexicaseSelect(*this, fit_set, POP_SIZE);
                 break;

            default:
                emp_assert(false && "INVALID SELECTION SCEHME", SELECTION);
                exit(1);
                break;
        }
        Update();
    }

    void Run() {
        for (size_t u = 0; u <= MAX_GENS; u++) {
            RunStep();
        }  
    }

    std::map<int, int> SimulateEcoEA(emp::vector<size_t> exclude, const emp::vector< std::function<double(const ORG &)> > & extra_funs,
                   emp::vector<emp::Resource> & resources, size_t t_size, size_t tourny_count=1, double frac = .0025, double max_bonus = 5, double cost = 0, bool use_base = true) {

        emp_assert(GetFitFun(), "Must define a base fitness function");
        emp_assert(GetSize() > 0);
        emp_assert(t_size > 0, t_size);

        std::map<int, int> selections;

        // Initialize map
        for (int i = 0; i < GetSize(); i++) {
            if (IsOccupied(i)) {
                selections[i] = 0;
            }
        }

        emp::vector<double> base_fitness(GetSize());

       // Collect all fitness info.
        for (size_t org_id = 0; org_id < GetSize(); org_id++) {
            if (!IsOccupied(org_id)) {
                continue;
            }

            if (emp::Has(exclude, org_id)) {
                base_fitness[org_id] = 0;
                continue; // Don't allow fitness gain from resources
            }

            if (use_base) {
                base_fitness[org_id] = CalcFitnessID(org_id);
            } else {
                base_fitness[org_id] = 0;
            }

            for (size_t ex_id = 0; ex_id < extra_funs.size(); ex_id++) {

                // resources[ex_id].Inc(resources[ex_id].GetInflow()/GetNumOrgs());
                double cur_fit = extra_funs[ex_id](GetOrg(org_id));
                cur_fit = emp::Pow(cur_fit, 2.0);            
                    cur_fit *= frac*(resources[ex_id].GetAmount()-cost);
                    if (cur_fit > 0) {
                        cur_fit -= cost;
                    } else {
                        cur_fit = 0;
                    }

                cur_fit = std::min(cur_fit, max_bonus);

                base_fitness[org_id] *= emp::Pow2(cur_fit);
                // resources[ex_id].Dec(std::abs(cur_fit));
            }
        }

        emp::vector<size_t> entries;
        for (size_t T = 0; T < tourny_count; T++) {
            entries.resize(0);
            for (size_t i=0; i<t_size; i++) entries.push_back( GetRandomOrgID() ); // Allows replacement!
            double best_fit = base_fitness[entries[0]];
            size_t best_id = entries[0];

            // Search for a higher fit org in the tournament.
            for (size_t i = 1; i < t_size; i++) {
                const double cur_fit = base_fitness[entries[i]];
                if (cur_fit > best_fit) {
                    best_fit = cur_fit;
                    best_id = entries[i];
                }
            }

            // Place the highest fitness into the next generation!
            selections[best_id]++;
        }

        return selections;
    }

    // This works for sharing too, since sharing is baked into the fitness function.
    std::map<int, int> SimulateTournament(emp::vector<size_t> exclude, size_t t_size, size_t tourny_count=1) {

        emp_assert(t_size > 0, "Cannot have a tournament with zero organisms.", t_size, GetNumOrgs());
        emp_assert(t_size <= GetNumOrgs(), "Tournament too big for world.", t_size, GetNumOrgs());
        emp_assert(tourny_count > 0);

        std::map<int, int> selections;

        // Initialize map
        for (int i = 0; i < GetSize(); i++) {
            if (IsOccupied(i)) {
                selections[i] = 0;
            }
        }

        emp::vector<size_t> entries;
        for (size_t T = 0; T < tourny_count; T++) {
            entries.resize(0);
            // Choose organisms for this tournament (with replacement!)
            for (size_t i=0; i < t_size; i++) entries.push_back( GetRandomOrgID() );


            double best_fit;
            if (emp::Has(exclude, entries[0])) {
                best_fit = 0;
            } else {
                best_fit = CalcFitnessID(entries[0]);
            }

            size_t best_id = entries[0];

            // Search for a higher fit org in the tournament.
            for (size_t i = 1; i < t_size; i++) {
                if (!emp::Has(exclude, entries[i])) {
                    const double cur_fit = CalcFitnessID(entries[i]);
                    if (cur_fit > best_fit) {
                        best_fit = cur_fit;
                        best_id = entries[i];
                    }
                }
            }

            // highest fitness wins!
            selections[best_id]++;
        }

        return selections;
    }


    EMP_CREATE_OPTIONAL_METHOD(TriggerOnLexicaseSelect, TriggerOnLexicaseSelect);
    std::map<int, int> SimulateLexicase(emp::vector<size_t> exclude, size_t repro_count=1, size_t max_funs = 0) {

        std::map<int, int> selections;

        // Initialize map
        for (int i = 0; i < GetSize(); i++) {
            if (IsOccupied(i)) {
                selections[i] = 0;
            }
        }


        std::map<typename ORG::genome_t, int> genotype_counts;
        emp::vector<emp::vector<size_t>> genotype_lists;

        // Find all orgs with same genotype - we can dramatically reduce
        // fitness evaluations this way.
        for (size_t org_id = 0; org_id < GetSize(); org_id++) {
            if (emp::Has(exclude, org_id)) {
                continue; // Excluded orgs will never be selected
            }

            if (IsOccupied(org_id)) {
                const auto gen = GetGenomeAt(org_id);
                if (emp::Has(genotype_counts, gen)) {
                    genotype_lists[genotype_counts[gen]].push_back(org_id);
                } else {
                    genotype_counts[gen] = genotype_lists.size();
                    genotype_lists.emplace_back(emp::vector<size_t>{org_id});
                }
            }
        }

        emp::vector<size_t> all_gens(genotype_lists.size()), cur_gens, next_gens;

        for (size_t i = 0; i < genotype_lists.size(); i++) {
            all_gens[i] = i;
        }

        if (!max_funs) max_funs = fit_set.size();
        // std::cout << "in lexicase" << std::endl;
        // Collect all fitness info. (@CAO: Technically only do this is cache is turned on?)
        emp::vector< emp::vector<double> > fitnesses(fit_set.size());
        for (size_t fit_id = 0; fit_id < fit_set.size(); ++fit_id) {
            fitnesses[fit_id].resize(genotype_counts.size());
            // std::cout << fit_id << std::endl;
            int id = 0;
            for (auto & gen : genotype_lists) {
                fitnesses[fit_id][id] = fit_set[fit_id](GetOrg(gen[0]));
                id++;
            }
        }

        // std::cout << to_string(fitnesses) << std::endl;
        // std::cout << "fitness calculated" << std::endl;
        // Go through a new ordering of fitness functions for each selections.
        // std::cout << "randdomized" << std::endl;
        for (size_t repro = 0; repro < repro_count; ++repro) {
            // std::cout << "repro: " << repro << std::endl;
            // Determine the current ordering of the functions.
            emp::vector<size_t> order;

            if (max_funs == fit_set.size()) {
                order = GetPermutation(GetRandom(), fit_set.size());
            } else {
                order.resize(max_funs); // We want to limit the total numebr of tests done.
                for (auto & x : order) x = GetRandom().GetUInt(fit_set.size());
            }
            // @CAO: We could have selected the order more efficiently!
            // std::cout << "reoreder" << std::endl;
            // Step through the functions in the proper order.
            cur_gens = all_gens;  // Start with all of the organisms.
            int depth = -1;
            for (size_t fit_id : order) {
                // std::cout << "fit_id: " << fit_id << std::endl;
                depth++;

                // std::cout << "about to index" << std::endl;
                // std::cout << to_string(fitnesses[fit_id]) << std::endl;
                // std::cout << cur_orgs[0] << std::endl;
                double max_fit = fitnesses[fit_id][cur_gens[0]];
                next_gens.push_back(cur_gens[0]);
                // std::cout << "Starting max: " << max_fit << to_string(cur_gens) << std::endl;
                for (size_t gen_id : cur_gens) {

                    const double cur_fit = fitnesses[fit_id][gen_id];
                    // std::cout << "gen_id: " << gen_id << "Fit: " << cur_fit << std::endl;
                    if (cur_fit > max_fit) {
                        max_fit = cur_fit;             // This is a the NEW maximum fitness for this function
                        next_gens.resize(0);           // Clear out orgs with former maximum fitness
                        next_gens.push_back(gen_id);   // Add this org as only one with new max fitness
                        // std::cout << "New max: " << max_fit << " " << cur_gens.size() << std::endl;
                    }
                    else if (cur_fit == max_fit) {
                        next_gens.push_back(gen_id);   // Same as cur max fitness -- save this org too.
                        // std::cout << "Adding: " << gen_id << std::endl;
                    }
                }
                // Make next_orgs into new cur_orgs; make cur_orgs allocated space for next_orgs.
                std::swap(cur_gens, next_gens);
                next_gens.resize(0);

                if (cur_gens.size() == 1) break;  // Stop if we're down to just one organism.
            }

            // Place a random survivor (all equal) into the next generation!
            emp_assert(cur_gens.size() > 0, cur_gens.size(), fit_set.size(), all_gens.size());
            size_t options = 0;
            for (size_t gen : cur_gens) {
                options += genotype_lists[gen].size();
            }
            size_t winner = GetRandom().GetUInt(options);
            int repro_id = -1;

            for (size_t gen : cur_gens) {
                if (winner < genotype_lists[gen].size()) {
                    repro_id = (int) genotype_lists[gen][winner];
                    break;
                }
                winner -= genotype_lists[gen].size();
            }
            emp_assert(repro_id != -1, repro_id, winner, options);

            selections[repro_id]++;

            // std::cout << depth << "abotu to calc used" <<std::endl;
            emp::vector<size_t> used = emp::Slice(order, 0, depth+1);
            // If the world has a OnLexicaseSelect method, call it
            // std::cout << depth << " " << to_string(used) << std::endl;
            TriggerOnLexicaseSelect(used, repro_id);
        }
        // std::cout << "Done with lex" << std::endl;
        return selections;
    } 

    emp::WeightedGraph CalcInteractionNetworks() {

        std::map<int, int> base_fitness;
        std::map<int, int> new_fitness;

        emp::WeightedGraph effects(GetSize());

        switch(SELECTION) {
            case (uint32_t)SELECTION_METHOD::TOURNAMENT :
                base_fitness = SimulateTournament({}, TOURNAMENT_SIZE, ECOLOGY_DATA_N_SIMULATIONS);
                break;

            case (uint32_t)SELECTION_METHOD::SHARING : // Sharing is handled in the setting of the fitness function
                base_fitness = SimulateTournament({}, TOURNAMENT_SIZE, ECOLOGY_DATA_N_SIMULATIONS);
                break;

            case (uint32_t)SELECTION_METHOD::ECOEA :
                base_fitness = SimulateEcoEA({}, fit_set, resources, TOURNAMENT_SIZE, ECOLOGY_DATA_N_SIMULATIONS, RESOURCE_SELECT_FRAC, RESOURCE_SELECT_MAX_BONUS,RESOURCE_SELECT_COST, true);
                break;

            case (uint32_t)SELECTION_METHOD::RANDOM :
                for (int i = 0; i < GetSize(); i++) {
                    if (IsOccupied(i)) {
                        base_fitness[i] = 1;
                    }
                }
                break;

            case (uint32_t)SELECTION_METHOD::LEXICASE :
                base_fitness = SimulateLexicase({}, fit_set, ECOLOGY_DATA_N_SIMULATIONS);
                break;

            default:
                emp_assert(false && "INVALID SELECTION SCEHME", SELECTION);
                exit(1);
                break;
        }

        for (int i = 0; i < GetSize(); i++) {
            if (!IsOccupied(i)) {
                continue;
            }

            switch(SELECTION) {
                case (uint32_t)SELECTION_METHOD::TOURNAMENT :
                    new_fitness = SimulateTournament({}, TOURNAMENT_SIZE, ECOLOGY_DATA_N_SIMULATIONS);
                    break;

                case (uint32_t)SELECTION_METHOD::SHARING : // Sharing is handled in the setting of the fitness function
                    new_fitness = SimulateTournament({}, TOURNAMENT_SIZE, ECOLOGY_DATA_N_SIMULATIONS);
                    break;

                case (uint32_t)SELECTION_METHOD::ECOEA :
                    new_fitness = SimulateEcoEA({}, fit_set, resources, TOURNAMENT_SIZE, ECOLOGY_DATA_N_SIMULATIONS, RESOURCE_SELECT_FRAC, RESOURCE_SELECT_MAX_BONUS,RESOURCE_SELECT_COST, true);
                    break;

                case (uint32_t)SELECTION_METHOD::RANDOM :
                    for (int i = 0; i < GetSize(); i++) {
                        if (IsOccupied(i)) {
                            new_fitness[i] = 1;
                        }
                    }
                    break;

                case (uint32_t)SELECTION_METHOD::LEXICASE :
                    new_fitness = SimulateLexicase({}, fit_set, ECOLOGY_DATA_N_SIMULATIONS);
                    break;

                default:
                    emp_assert(false && "INVALID SELECTION SCEHME", SELECTION);
                    exit(1);
                    break;
            }

            for (auto res : new_fitness) {
                effects.AddEdge(i, res.first, base_fitness[res.first] - new_fitness[res.first]/(double)POP_SIZE);
            }
            
        }

        return effects;
    }

    void ComputeNetworks() {
        emp::WeightedGraph i_network = CalcInteractionNetworks();

        for (emp::vector<double> & weights : i_network.GetWeights()) {
            for (double w : weights) {
                interaction_strengths.Add(w);
            }
        }

        for (size_t i = 0; i < i_network.GetSize(); i++) {
            node_out_degrees.Add(i_network.GetDegree(i));
            node_in_degrees.Add(i_network.GetInDegree(i));
        }


    }

};