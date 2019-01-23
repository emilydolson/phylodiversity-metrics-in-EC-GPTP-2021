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

#include "interaction_networks.h"

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
  VALUE(RESOURCE_SELECT_NICHE_WIDTH, double, 0.0, "Score necessary to count as attempting to use a resource"),
  VALUE(LEXICASE_EPSILON, double, 0.0, "To use epsilon-lexicase, set this to an epsilon value greater than 0"),

  GROUP(OPEN_ENDED_EVOLUTION, "Settings related to tracking MODES metrics"),
  VALUE(MODES_RESOLUTION, int, 1, "How often should MODES metrics be calculated?"),
  VALUE(FILTER_LENGTH, int, 1000, "How many generations should we use for the persistence filter?"),

  GROUP(DATA_RESOLUTION, "How often should data get printed?"),
  VALUE(FAST_DATA_RES, int, 10, "How often should easy to calculate metrics (e.g. fitness) be calculated?"),
  VALUE(ECOLOGY_DATA_RES, int, 100, "How often should ecological interactions (expensive) be calculated?")
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

enum class PROBLEM_TYPE { NK=0, PROGRAM_SYNTHESIS=1, REAL_VALUE=2, SORTING_NETWORK=3, LOGIC_9=4 };

template <typename ORG>
class EcologyWorld : public emp::World<ORG> {

    using phen_t = emp::vector<double>;
    using fit_map_t = std::map<phen_t, double>;

    using typename emp::World<ORG>::fun_calc_fitness_t;
    using typename emp::World<ORG>::genome_t;
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
    using emp::World<ORG>::AddSystematics;
    using emp::World<ORG>::OnUpdate;
    using emp::World<ORG>::OnPlacement;
    using emp::World<ORG>::OnBeforePlacement;
    using emp::World<ORG>::SetMutFun;
    using emp::World<ORG>::Inject;
    using emp::World<ORG>::SetFitFun;
    using emp::World<ORG>::SetSharedFitFun;
    using emp::World<ORG>::SetupFile;
    using emp::World<ORG>::GetFitFun;
    using emp::World<ORG>::GetNumOrgs;

    enum class SELECTION_METHOD { TOURNAMENT=0, SHARING=1, LEXICASE=2, ECOEA=3, RANDOM=4 };

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
    int FAST_DATA_RES;
    int ECOLOGY_DATA_RES;

    double SHARING_ALPHA;
    double SHARING_THRESHOLD;
    double RESOURCE_SELECT_RES_AMOUNT;
    double RESOURCE_SELECT_RES_INFLOW;
    double RESOURCE_SELECT_RES_OUTFLOW;
    double RESOURCE_SELECT_FRAC;
    double RESOURCE_SELECT_MAX_BONUS;
    double RESOURCE_SELECT_COST;
    double RESOURCE_SELECT_NICHE_WIDTH;
    double LEXICASE_EPSILON;

    emp::NKLandscape landscape;

    // Function to determine actual performance on problem
    fun_calc_fitness_t evaluation_fun;
    std::function<fit_map_t(emp::vector<phen_t> &, all_attrs)> evaluate_competition_fun;
    all_attrs attrs = DEFAULT;

    emp::vector<phen_t> phenotypes;

    // Set of fitness functions for lexicase and Eco-EA
    emp::vector<fun_calc_fitness_t> fit_set;
    emp::vector<emp::Resource> resources;

    emp::Ptr<emp::OEETracker<emp::Systematics<ORG, ORG>, ORG>> oee;

    double interaction_signif_threshold = 0;

    // Data tracking
    emp::DataNode<double, emp::data::Histogram, emp::data::Stats> pos_interaction_strengths;
    emp::DataNode<double, emp::data::Histogram, emp::data::Stats> neg_interaction_strengths;
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
        RESOURCE_SELECT_NICHE_WIDTH = config.RESOURCE_SELECT_NICHE_WIDTH();
        LEXICASE_EPSILON = config.LEXICASE_EPSILON();

        attrs = emp::tools::Merge(SigmaShare(SHARING_THRESHOLD), Alpha(SHARING_ALPHA),
                                  Cost(RESOURCE_SELECT_COST), Cf(RESOURCE_SELECT_FRAC),
                                  NicheWidth(RESOURCE_SELECT_NICHE_WIDTH), MaxScore(1),
                                  ResourceInflow(RESOURCE_SELECT_RES_INFLOW), 
                                  ResourceOutflow(RESOURCE_SELECT_RES_OUTFLOW), 
                                  MaxBonus(RESOURCE_SELECT_MAX_BONUS), 
                                  TournamentSize(TOURNAMENT_SIZE), 
                                  Epsilon(LEXICASE_EPSILON));

        emp::Ptr<emp::Systematics<ORG, ORG> > sys;
        sys.New([](const ORG & o){return o;});
        oee.New(sys, [this](ORG & org){return org;}, [](const ORG & org){
            // TODO: Actually calculate complexity
            return 1;
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
            perf_node->Reset();
            for (emp::Ptr<ORG> org : pop) {
                perf_node->Add(evaluation_fun(*org));
            }
        });

        OnUpdate([this](size_t ud){if(ud % FAST_DATA_RES == 0) {evaluate_performance_sig.Trigger();}});
        OnUpdate([this](size_t ud){if(ud % ECOLOGY_DATA_RES == 0) {ComputeNetworks();}});

        emp::DataFile & performance_file = SetupFile("performance.csv");
        performance_file.AddVar(update, "update", "Update");
        performance_file.AddStats(*GetDataNode("performance"), "performance", "How well are programs actually solving the problem?");
        performance_file.PrintHeaderKeys();
        performance_file.SetTimingRepeat(FAST_DATA_RES);


        emp::DataFile & oee_file = SetupFile("oee.csv");
        oee_file.AddVar(update, "generation", "Generation");
        oee_file.AddCurrent(*oee->GetDataNode("change"), "change", "change potential");
        oee_file.AddCurrent(*oee->GetDataNode("novelty"), "novelty", "novelty potential");
        oee_file.AddCurrent(*oee->GetDataNode("diversity"), "ecology", "ecology potential");
        oee_file.AddCurrent(*oee->GetDataNode("complexity"), "complexity", "complexity potential");
        oee_file.PrintHeaderKeys();
        oee_file.SetTimingRepeat(MODES_RESOLUTION);

        pos_interaction_strengths.SetupBins(0.0, 1.01, 20);
        neg_interaction_strengths.SetupBins(-1.0, 0.0, 20);
        node_out_degrees.SetupBins(0, POP_SIZE+1, 20);
        node_in_degrees.SetupBins(0, POP_SIZE+1, 20);

        emp::DataFile & ecology_file = SetupFile("ecology.csv");
        ecology_file.AddVar(update, "generation", "Generation");
        ecology_file.AddStats(pos_interaction_strengths, "pos_interaction_strength", "positive interaction strength");
        ecology_file.AddAllHistBins(pos_interaction_strengths, "pos_interaction_strength", "positive interaction strength");
        ecology_file.AddStats(neg_interaction_strengths, "neg_interaction_strength", "negative interaction strength");
        ecology_file.AddAllHistBins(neg_interaction_strengths, "neg_interaction_strength", "negative interaction strength");
        ecology_file.AddStats(node_out_degrees, "interaction_node_out_degree", "interaction node out degree"); 
        ecology_file.AddAllHistBins(node_out_degrees, "interaction_node_out_degree", "interaction node out degree"); 
        ecology_file.AddStats(node_in_degrees, "interaction_node_in_degree", "interaction node in degree");
        ecology_file.AddAllHistBins(node_in_degrees, "interaction_node_in_degree", "interaction node in degree");
        ecology_file.PrintHeaderKeys();
        ecology_file.SetTimingRepeat(ECOLOGY_DATA_RES);

        switch (PROBLEM){
            case (uint32_t) PROBLEM_TYPE::NK:
                if constexpr (std::is_same<ORG, emp::BitVector>::value) {
                    // We need this if so we don't need to define SetupNK for every
                    // possible ORG type
                    SetupNK();
                } else {
                    exit(1);
                }
                break;

            case (uint32_t) PROBLEM_TYPE::PROGRAM_SYNTHESIS:
                // TODO
                break;            

            case (uint32_t) PROBLEM_TYPE::REAL_VALUE:
                // TODO
                break;            

            case (uint32_t) PROBLEM_TYPE::SORTING_NETWORK:
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

            case ( (uint32_t) SELECTION_METHOD::TOURNAMENT):
                evaluate_competition_fun = do_tournament<phen_t>;
                break;

            case ( (uint32_t) SELECTION_METHOD::SHARING):
                evaluate_competition_fun = do_sharing<phen_t>;
                break;

            case ( (uint32_t) SELECTION_METHOD::LEXICASE):
                evaluate_competition_fun = do_lexicase<phen_t>;
                break;

            case ( (uint32_t) SELECTION_METHOD::ECOEA):

                evaluate_competition_fun = do_eco_ea<phen_t>;

                resources.resize(0);
                for (size_t i=0; i<N_TEST_CASES; i++) {
                    resources.push_back(emp::Resource(RESOURCE_SELECT_RES_INFLOW, RESOURCE_SELECT_RES_INFLOW, .01));
                }
                break;

            default:
                break;

        }

        SetAutoMutate();
        if (SELECTION == (uint32_t)SELECTION_METHOD::TOURNAMENT) {
            SetCache();
        }

        phenotypes.resize(POP_SIZE);

        OnUpdate([this](int ud){phenotypes.clear(); phenotypes.resize(POP_SIZE);});
        OnBeforePlacement([this](ORG & o, int pos){
            for (auto fun : fit_set) {
                phenotypes[pos].push_back(fun(o));
            }
        }); 

        InitPop();
    }


    void InitPop();

    void SetupNK();

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
                emp::ResourceSelect(*this, fit_set, resources, TOURNAMENT_SIZE, POP_SIZE, RESOURCE_SELECT_FRAC, RESOURCE_SELECT_MAX_BONUS,RESOURCE_SELECT_COST, true, RESOURCE_SELECT_NICHE_WIDTH);
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

    void ComputeNetworks() {
        emp::WeightedGraph i_network = CalcCompetition(phenotypes, evaluate_competition_fun, attrs);
        neg_interaction_strengths.Reset();
        pos_interaction_strengths.Reset();
        node_out_degrees.Reset();
        node_in_degrees.Reset();

        for (emp::vector<double> & weights : i_network.GetWeights()) {
            for (double w : weights) {
                if (w < 0){
                    neg_interaction_strengths.Add(w);
                } else if (w > 0) {
                    pos_interaction_strengths.Add(w);
                }
            }
        }

        for (size_t i = 0; i < i_network.GetSize(); i++) {
            node_out_degrees.Add(i_network.GetDegree(i));
            node_in_degrees.Add(i_network.GetInDegree(i));
        }


    }

};

template <>
void EcologyWorld<emp::BitVector>::InitPop() {
    // Build a random initial population
    for (uint32_t i = 0; i < POP_SIZE; i++) {
        emp::BitVector next_org(N);
        for (uint32_t j = 0; j < N; j++) next_org[j] = random_ptr->P(0.5);
        Inject(next_org);
    }
}

template <>
void EcologyWorld<emp::vector<double>>::InitPop() {
    // Build a random initial population
    for (uint32_t i = 0; i < POP_SIZE; i++) {
        emp::vector<double> next_org(N);
        for (uint32_t j = 0; j < N; j++) next_org[j] = random_ptr->GetDouble(1);
        Inject(next_org);
    }
}


template <>
void EcologyWorld<emp::BitVector>::SetupNK() {
        // Create landscape based on correct random seed
        landscape = emp::NKLandscape(N, K, *random_ptr);

        // Set-up fitness function (queries NK Landscape)
        evaluation_fun =
            [this](emp::BitVector & org){ return landscape.GetFitness(org); };

        if (SELECTION == (uint32_t) SELECTION_METHOD::SHARING) {
            SetSharedFitFun(evaluation_fun, [](emp::BitVector & org1, emp::BitVector & org2){return calc_hamming_distance(org1, org2);}, SHARING_THRESHOLD, SHARING_ALPHA);
        } else {
            std::cout << "Setting fit fun" << std::endl;
            SetFitFun(evaluation_fun);
        }

        // Make a fitness function for each gene (set of sites that determine a fitness contribution)
        // in the bitstring
 
        for (size_t gene = 0; gene < N; gene++) {
            fit_set.push_back([this, gene](emp::BitVector & org){return landscape.GetFitness(gene, org);});
        }

        // Setup the mutation function.
        std::function<size_t(emp::BitVector &, emp::Random &)> mut_fun =
            [this](emp::BitVector & org, emp::Random & random) {
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