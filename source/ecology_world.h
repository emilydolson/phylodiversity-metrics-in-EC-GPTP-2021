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
#include "tools/Random.h"
#include "tools/sequence_utils.h"
#include "Evolve/OEE.h"
#include "tools/vector_utils.h"
#include "tools/Graph.h"

#include "org_types.h"
#include "TaskSet.h"
#include "TestcaseSet.h"
#include "BitSorterMutators.h"

#include "interaction_networks.h"

#include "cec2013.h"

EMP_BUILD_CONFIG( EcologyConfig,
  GROUP(MAIN, "Global settings"),
  VALUE(MODE, int, 0, "0 = run experiment, 1 = Just run ecology stats so we can time it."),
  VALUE(SEED, int, 0, "Random number seed (0 for based on time)"),
  VALUE(START_POP_SIZE, uint32_t, 1, "Number of organisms to seed population with."),
  VALUE(POP_SIZE, uint32_t, 1000, "Number of organisms in the popoulation."),
  VALUE(MAX_GENS, uint32_t, 2000, "How many generations should we process?"),
  VALUE(MUT_RATE, double, .005, "Probability of each site being mutated. For real-valued problems, the standard deviation of of the distribution from which mutations are pulled. For program, the probability of instruction mutation to a different one."),
  VALUE(PROBLEM, uint32_t, 0, "Which problem to use? 0 = NK, 1 = Program synthesis, 2 = Real-valued, 3 = Sorting network, 4 = Logic-9"),

  GROUP(NK, "Settings for NK landscape"),
  VALUE(K, uint32_t, 10, "Level of epistasis in the NK model"),
  VALUE(N, uint32_t, 200, "Number of bits in each organisms (must be > K)"),

  GROUP(REAL_VALUED, "Settings for real-valued optimzation problems"),
  VALUE(FUNCTION_NUMBER, uint32_t, 0, "Problem to use"),
  VALUE(DIMS, uint32_t, 200, "Number of dimensions in orgs"),

  GROUP(SORTING_NETWORKS, "Setting for sorting network problems"),
  VALUE(NUM_BITS, size_t, 16, "Length of input to sorting network."),
  VALUE(MAX_NETWORK_SIZE, size_t, 128, "Maximum size of a sorting network."),
  VALUE(MIN_NETWORK_SIZE, size_t, 1, "Minimum size of a sorting network."),
  VALUE(PER_INDEX_SUB, double, 0.001, "."),
  VALUE(PER_PAIR_DUP, double, 0.0005, "."),
  VALUE(PER_PAIR_INS, double, 0.0005, "."),
  VALUE(PER_PAIR_DEL, double, 0.001, "."),
  VALUE(PER_PAIR_SWAP, double, 0.001, "."),

  GROUP(PROGRAM_SYNTHESIS, "Settings for evolving code"),
  VALUE(GENOME_SIZE, uint32_t, 100, "Starting length of genome"),
  VALUE(MAX_GENOME_SIZE, uint32_t, 1000, "Maximum length of genome"),
  VALUE(EVAL_TIME, uint32_t, 100, "How long to run program for"),
  VALUE(INS_MUT_RATE, double, 0.005, "Probability of insertion mutation per instruction"),  
  VALUE(DEL_MUT_RATE, double, 0.005, "Probability of deletion mutation per instruction"),  
  VALUE(ARG_SUB_RATE, double, 0.005, "Probability argument substituion per instruction"),  

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
emp::vector<bool> Skeletonize(bit_t & org, std::function<double(bit_t&)> fun_calc_fitness) {
    emp_assert(org.size() > 0, "Empty org passed to skeletonize");

    emp::vector<bool> skeleton;
    double fitness = fun_calc_fitness(org);
    bit_t test_org = bit_t(org);

    for (size_t i = 0; i < org.size(); i++) {
        test_org[i] = !test_org[i]; // For bit strings we just flip the bit
        double new_fitness = fun_calc_fitness(test_org);
        if (new_fitness < fitness) {
            skeleton.push_back(true);
        } else {
            skeleton.push_back(false);
        }
        test_org[i] = (int)org[i];
    }

    return skeleton;
}

const emp::vector<size_t> PROBLEM_MAP = {4,5,6,7,10,11,12,13};
const emp::vector<std::string> PROBLEM_DESC = {"F4 (2D)","F5 (2D)","F6 (2D)","F7 (2D)","F8 (2D)","F9 (2D)","F10 (2D)","F11 (2D)"};

enum class PROBLEM_TYPE { NK=0, PROGRAM_SYNTHESIS=1, REAL_VALUE=2, SORTING_NETWORK=3, LOGIC_9=4 };

template <typename ORG>
class EcologyWorld : public emp::World<ORG> {

    using phen_t = emp::vector<double>;
    using fit_map_t = emp::vector<double>;
    using phen_taxon_t = emp::Taxon<phen_t, emp::datastruct::mut_landscape_info<phen_t>>;

    // Logic-9 constants
    static constexpr double MIN_POSSIBLE_SCORE = -32767;
    static constexpr size_t MAX_LOGIC_TASK_NUM_INPUTS = 2;
    static constexpr uint32_t MIN_LOGIC_TASK_INPUT = 0;
    static constexpr uint32_t MAX_LOGIC_TASK_INPUT = 1000000000;

    using task_io_t = uint32_t;
    using taskset_t = TaskSet<std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS>, task_io_t>;
    using test_case_t = std::pair<emp::vector<int>, double>;

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
    uint32_t START_POP_SIZE;
    uint32_t EVAL_TIME;
    uint32_t MAX_GENS;
    uint32_t SELECTION;
    uint32_t CHANGE_TYPE;
    uint32_t CHANGE_RATE;
    uint32_t PROBLEM;
    uint32_t N_TEST_CASES;
    uint32_t GENOME_SIZE;
    uint32_t MAX_GENOME_SIZE;
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

    double INS_MUT_RATE;
    double DEL_MUT_RATE;
    double ARG_SUB_RATE;
    double PER_INDEX_SUB;
    double PER_PAIR_DUP;
    double PER_PAIR_INS;
    double PER_PAIR_DEL;
    double PER_PAIR_SWAP;

    size_t NUM_BITS;
    size_t MAX_NETWORK_SIZE;
    size_t MIN_NETWORK_SIZE;

    std::string TESTSET_FILE;
    std::string TRAINSET_FILE;

    // Problem-specific member variables
    emp::NKLandscape landscape; // NK
    taskset_t task_set; // Logic-9
    std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS> task_inputs; // Logic-9
    size_t input_load_id; // Logic-9
    gp_t::inst_lib_t inst_set; // AvidaGP problems
    TestcaseSet<int, double> training_set; // Program synthesis
    TestcaseSet<int, double> testing_set; // Program synthesis
    BitSorterMutator sorter_mutator;

    // Real value
    emp::Ptr<CEC2013> optim_function;
    emp::vector<emp::vector<double>> key_points;
    emp::vector<double> mid_point;
    emp::vector<double> ubounds;
    emp::vector<double> lbounds;
    emp::vector<double> mut_std;


    std::map<ORG, phen_t> phen_map;

    // Function to determine actual performance on problem
    fun_calc_fitness_t evaluation_fun;
    std::function<fit_map_t(emp::vector<phen_t> &, all_attrs)> evaluate_competition_fun;
    all_attrs attrs = DEFAULT;

    // Set of fitness functions for lexicase and Eco-EA
    emp::vector<fun_calc_fitness_t> fit_set;
    emp::vector<emp::Resource> resources;

    emp::Ptr<emp::OEETracker<emp::Systematics<ORG, ORG>, emp::vector<typename ORG::gene_t> > > oee;


    // Data tracking
    emp::DataNode<double, emp::data::Histogram, emp::data::Stats> pos_interaction_strengths;
    emp::DataNode<double, emp::data::Histogram, emp::data::Stats> neg_interaction_strengths;
    emp::DataNode<int, emp::data::Histogram, emp::data::Stats> node_in_degrees;
    emp::DataNode<int, emp::data::Histogram, emp::data::Stats> node_out_degrees;
    emp::Ptr<emp::Systematics<ORG, ORG> > sys;
    emp::Ptr<emp::Systematics<ORG, phen_t, emp::datastruct::mut_landscape_info<phen_t>> > phen_sys;

    emp::Ptr<ORG> best_curr;
    double best_curr_fit = 0;

    public:

    EcologyWorld() {;}
    EcologyWorld(emp::Random & rnd) : emp::World<ORG>(rnd) {;}
    ~EcologyWorld(){if (oee){oee.Delete();}}

    void Setup(EcologyConfig & config) {
        N = config.N();
        K = config.K();
        POP_SIZE = config.POP_SIZE();
        START_POP_SIZE = config.START_POP_SIZE();
        EVAL_TIME = config.EVAL_TIME();
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
        GENOME_SIZE = config.GENOME_SIZE();
        MAX_GENOME_SIZE = config.MAX_GENOME_SIZE();

        RESOURCE_SELECT_RES_AMOUNT = config.RESOURCE_SELECT_RES_AMOUNT();
        RESOURCE_SELECT_RES_INFLOW = config.RESOURCE_SELECT_RES_INFLOW();
        RESOURCE_SELECT_RES_OUTFLOW = config.RESOURCE_SELECT_RES_OUTFLOW();
        RESOURCE_SELECT_FRAC = config.RESOURCE_SELECT_FRAC();
        RESOURCE_SELECT_MAX_BONUS = config.RESOURCE_SELECT_MAX_BONUS();
        RESOURCE_SELECT_COST = config.RESOURCE_SELECT_COST();
        RESOURCE_SELECT_NICHE_WIDTH = config.RESOURCE_SELECT_NICHE_WIDTH();
        LEXICASE_EPSILON = config.LEXICASE_EPSILON();

        INS_MUT_RATE = config.INS_MUT_RATE();
        DEL_MUT_RATE = config.DEL_MUT_RATE();
        ARG_SUB_RATE = config.ARG_SUB_RATE();
        PER_INDEX_SUB = config.PER_INDEX_SUB();
        PER_PAIR_DUP = config.PER_PAIR_DUP();
        PER_PAIR_INS = config.PER_PAIR_INS();
        PER_PAIR_DEL = config.PER_PAIR_DEL();
        PER_PAIR_SWAP = config.PER_PAIR_SWAP();
  
        MAX_NETWORK_SIZE = config.MAX_NETWORK_SIZE();
        MIN_NETWORK_SIZE = config.MIN_NETWORK_SIZE();
        NUM_BITS = config.NUM_BITS();

        TESTSET_FILE = config.TESTSET_FILE();
        TRAINSET_FILE = config.TRAINSET_FILE();

        attrs = emp::tools::Merge(SigmaShare(SHARING_THRESHOLD), Alpha(SHARING_ALPHA),
                                  Cost(RESOURCE_SELECT_COST), Cf(RESOURCE_SELECT_FRAC),
                                  NicheWidth(RESOURCE_SELECT_NICHE_WIDTH), MaxScore(1),
                                  ResourceInflow(RESOURCE_SELECT_RES_INFLOW), 
                                  ResourceOutflow(RESOURCE_SELECT_RES_OUTFLOW), 
                                  MaxBonus(RESOURCE_SELECT_MAX_BONUS), 
                                  TournamentSize(TOURNAMENT_SIZE), 
                                  Epsilon(LEXICASE_EPSILON));

        sys.New([](const ORG & o){return o;});
        oee.New(sys, [this](ORG & org){
                if constexpr (std::is_same<ORG, bit_t>::value) {
                    return Skeletonize(org, GetFitFun());                    
                } else if constexpr (std::is_same<ORG, rv_t>::value) {
                    return org;
                } else if constexpr (std::is_same<ORG, gp_t>::value) {
                    gp_t::Instruction null_inst(inst_set.GetID("Null"));
                    return emp::Skeletonize(org, null_inst, GetFitFun());
                } else if constexpr (std::is_same<ORG, sorting_t>::value) {
                    emp::vector<sorting_t::gene_t> skel;
                    double fit = fun_calc_fitness(org);
                    for (int i = org.GetSize() - 1; i >= 0; i--) {
                        sorting_t test_org = org;
                        test_org.RemoveCompare(i);
                        if (fun_calc_fitness(test_org) < fit) {
                            skel.push_back(org.GetBits(i));
                        }
                    }
                    return skel;
                } else {
                    return org;
                }
            
            }, [](const emp::vector<typename ORG::gene_t> & org){
                if constexpr (std::is_same<ORG, bit_t>::value) {
                    return emp::Count(org, true);                    
                } else if constexpr (std::is_same<ORG, rv_t>::value) {
                    return -1;
                } else {
                    return org.size();
                }
        });
        oee->SetResolution(MODES_RESOLUTION);
        oee->SetGenerationInterval(FILTER_LENGTH);
        AddSystematics(sys);
        OnUpdate([this](int ud){oee->Update(ud);});

        phen_sys.New([this](ORG & o){            
            phen_t phen;
            for (auto fun : fit_set) {
                phen.push_back(fun(o));
            }
            return phen;});

        AddSystematics(phen_sys);
        std::function<void(emp::Ptr<phen_taxon_t>, ORG&)> record_taxon_data = [](emp::Ptr<phen_taxon_t> tax, ORG & o){
            tax->GetData().RecordFitness(o.fitness);
            tax->GetData().RecordPhenotype(o.phenotype);
        };
        phen_sys->OnNew(record_taxon_data);

        SetupFitnessFile().SetTimingRepeat(10);
        SetupSystematicsFile().SetTimingRepeat(10);
        SetupPopulationFile().SetTimingRepeat(10);
        SetPopStruct_Mixed(true);
        SetSynchronousSystematics(true);

        AddDataNode("performance");
        evaluate_performance_sig.AddAction([this](){
            auto perf_node = GetDataNode("performance");
            perf_node->Reset();
            best_curr = pop[0];
            best_curr_fit = 0.0;

            for (emp::Ptr<ORG> org : pop) {
                double fit = evaluation_fun(*org);
                perf_node->Add(fit);
                if (fit > best_curr_fit) {
                    best_curr = org;
                    best_curr_fit = fit;
                }
            }

            if constexpr (std::is_same<ORG, rv_t>::value) {
                std::cout << emp::to_string(*best_curr) << std::endl;
            } else {
                std::cout << *best_curr << std::endl;
            }
        });

        emp::DataFile & performance_file = SetupFile("performance.csv");
        performance_file.AddPreFun([this](){evaluate_performance_sig.Trigger();});
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
        ecology_file.AddPreFun([this](){ComputeNetworks();});
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

        emp::DataFile & species_ecology_file = SetupFile("species_ecology.csv");
        species_ecology_file.AddPreFun([this](){ComputeNetworks(true);});
        species_ecology_file.AddVar(update, "generation", "Generation");
        species_ecology_file.AddStats(pos_interaction_strengths, "pos_interaction_strength", "positive interaction strength");
        species_ecology_file.AddAllHistBins(pos_interaction_strengths, "pos_interaction_strength", "positive interaction strength");
        species_ecology_file.AddStats(neg_interaction_strengths, "neg_interaction_strength", "negative interaction strength");
        species_ecology_file.AddAllHistBins(neg_interaction_strengths, "neg_interaction_strength", "negative interaction strength");
        species_ecology_file.AddStats(node_out_degrees, "interaction_node_out_degree", "interaction node out degree"); 
        species_ecology_file.AddAllHistBins(node_out_degrees, "interaction_node_out_degree", "interaction node out degree"); 
        species_ecology_file.AddStats(node_in_degrees, "interaction_node_in_degree", "interaction node in degree");
        species_ecology_file.AddAllHistBins(node_in_degrees, "interaction_node_in_degree", "interaction node in degree");
        species_ecology_file.PrintHeaderKeys();
        species_ecology_file.SetTimingRepeat(ECOLOGY_DATA_RES);

        emp::DataFile & phylodiversity_file = SetupFile("phylodiversity.csv");
        sys->AddEvolutionaryDistinctivenessDataNode();
        sys->AddPairwiseDistanceDataNode();
        sys->AddPhylogeneticDiversityDataNode();

        phen_sys->AddEvolutionaryDistinctivenessDataNode();
        phen_sys->AddPairwiseDistanceDataNode();
        phen_sys->AddPhylogeneticDiversityDataNode();

        phylodiversity_file.AddVar(update, "generation", "Generation");
        phylodiversity_file.AddStats(*sys->GetDataNode("evolutionary_distinctiveness") , "genotype_evolutionary_distinctiveness", "evolutionary distinctiveness for a single update", true, true);
        phylodiversity_file.AddStats(*sys->GetDataNode("pairwise_distance"), "genotype_pairwise_distance", "pairwise distance for a single update", true, true);
        phylodiversity_file.AddCurrent(*sys->GetDataNode("phylogenetic_diversity"), "genotype_current_phylogenetic_diversity", "current phylogenetic_diversity", true, true);
        phylodiversity_file.AddStats(*phen_sys->GetDataNode("evolutionary_distinctiveness") , "phenotype_evolutionary_distinctiveness", "evolutionary distinctiveness for a single update", true, true);
        phylodiversity_file.AddStats(*phen_sys->GetDataNode("pairwise_distance"), "phenotype_pairwise_distance", "pairwise distance for a single update", true, true);
        phylodiversity_file.AddCurrent(*phen_sys->GetDataNode("phylogenetic_diversity"), "phenotype_current_phylogenetic_diversity", "current phylogenetic_diversity", true, true);
        phylodiversity_file.PrintHeaderKeys();
        phylodiversity_file.SetTimingRepeat(ECOLOGY_DATA_RES);

        emp::DataFile & dominant_file = SetupFile("dominant.csv");
        dominant_file.AddVar(update, "generation", "Generation");
        dominant_file.AddVar(best_curr_fit, "dominant_fitness", "Fitness of the best org in the population");
        dominant_file.AddFun<int>([this](){return best_curr->size();}, "dominant_size", "Size of dominant");
        dominant_file.AddFun<int>([this](){return emp::LineageLength(phen_sys->GetTaxonAt(best_curr->pos));}, "dominant_phenotypic_volatilty", "Number of phenotypes on dominant's lineage");
        dominant_file.AddFun<int>([this](){return emp::LineageLength(sys->GetTaxonAt(best_curr->pos));}, "dominant_lin_length", "Number of genotypes on dominant's lineage");
        dominant_file.AddFun<int>([this](){return emp::CountDeleteriousSteps(phen_sys->GetTaxonAt(best_curr->pos));}, "dominant_deleterious_steps", "Number of deleterious steps on dominant's lineage");
        dominant_file.AddFun<int>([this](){return emp::CountUniquePhenotypes(phen_sys->GetTaxonAt(best_curr->pos));}, "dominant_unique_phenotypes", "Number of unique phenotypes on dominant's lineage");
        dominant_file.PrintHeaderKeys();
        dominant_file.SetTimingRepeat(FAST_DATA_RES);

        emp::DataFile & lineage_file = SetupFile("lineage.csv");
        phen_sys->AddDeleteriousStepDataNode();
        phen_sys->AddVolatilityDataNode();
        phen_sys->AddUniqueTaxaDataNode();

        lineage_file.AddStats(*phen_sys->GetDataNode("deleterious_steps"), "deleterious_steps", "counts of deleterious steps along each lineage", true, true);
        lineage_file.AddStats(*phen_sys->GetDataNode("volatility"), "taxon_volatility", "counts of changes in taxon along each lineage", true, true);
        lineage_file.AddStats(*phen_sys->GetDataNode("unique_taxa"), "unique_taxa", "counts of unique taxa along each lineage", true, true);
        lineage_file.PrintHeaderKeys();
        lineage_file.SetTimingRepeat(FAST_DATA_RES);

        switch (PROBLEM){
            case (uint32_t) PROBLEM_TYPE::NK:
                if constexpr (std::is_same<ORG, bit_t>::value) {
                    // We need this if so we don't need to define SetupNK for every
                    // possible ORG type
                    SetupNK();
                } else {
                    exit(1);
                }
                break;

            case (uint32_t) PROBLEM_TYPE::PROGRAM_SYNTHESIS:
                if constexpr (std::is_same<ORG, gp_t>::value) {
                    SetupProgramSynthesis();
                } else {
                    exit(1);
                }
                break;            

            case (uint32_t) PROBLEM_TYPE::REAL_VALUE:
                // TODO
                break;            

            case (uint32_t) PROBLEM_TYPE::SORTING_NETWORK:
                if constexpr (std::is_same<ORG, sorting_t>::value) {
                    SetupSortingNetworks();
                } else {
                    exit(1);
                }
                break;            

            case (uint32_t) PROBLEM_TYPE::LOGIC_9:
                if constexpr (std::is_same<ORG, gp_t>::value) {
                    SetupLogic9();
                } else {
                    exit(1);
                }
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

        OnBeforePlacement([this](ORG & o, int pos){
            o.phenotype.resize(0);
            for (auto fun : fit_set) {
                o.phenotype.push_back(fun(o));
            }
            o.pos = pos;
            o.fitness = emp::Sum(o.phenotype);
        }); 

        InitPop();
        best_curr = pop[0];
    }


    void InitPop();
    void SetupNK();
    void SetupProgramSynthesis();
    void SetupLogic9();
    void SetupRealValue();
    void SetupSortingNetworks();
    double RunTestcase(gp_t & org, test_case_t testcase);

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
                emp::OptimizedLexicaseSelect(*this, fit_set, POP_SIZE);
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

    void TimeCompetition() {
        START_POP_SIZE = POP_SIZE;
        pop.clear();
        InitPop();
        emp::vector<phen_t> phenotypes;
        for (emp::Ptr<ORG> org : pop) {
            phenotypes.push_back(org->phenotype);
        }
        CalcCompetition(phenotypes, evaluate_competition_fun, attrs);        
    }

    void ComputeNetworks(bool species_level = false) {
        emp::vector<phen_t> phenotypes;
        for (emp::Ptr<ORG> org : pop) {
            phenotypes.push_back(org->phenotype);
        }

        if (species_level) {
            phenotypes = emp::RemoveDuplicates(phenotypes);
        }

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

    // For logic-9
    void ResetTasks() {
        input_load_id = 0;
        task_inputs[0] = random_ptr->GetUInt(MIN_LOGIC_TASK_INPUT, MAX_LOGIC_TASK_INPUT);
        task_inputs[1] = random_ptr->GetUInt(MIN_LOGIC_TASK_INPUT, MAX_LOGIC_TASK_INPUT);
        task_set.SetInputs(task_inputs);
        while (task_set.IsCollision()) {
            task_inputs[0] = random_ptr->GetUInt(MIN_LOGIC_TASK_INPUT, MAX_LOGIC_TASK_INPUT);
            task_inputs[1] = random_ptr->GetUInt(MIN_LOGIC_TASK_INPUT, MAX_LOGIC_TASK_INPUT);
            task_set.SetInputs(task_inputs);
        }
    }

};

template <>
void EcologyWorld<sorting_t>::InitPop() {
    sorter_mutator.MAX_NETWORK_SIZE = MAX_NETWORK_SIZE;
    sorter_mutator.MIN_NETWORK_SIZE = MIN_NETWORK_SIZE;
    sorter_mutator.SORT_SEQ_SIZE = NUM_BITS;
    sorter_mutator.PER_INDEX_SUB = PER_INDEX_SUB;
    sorter_mutator.PER_PAIR_DUP = PER_PAIR_DUP;
    sorter_mutator.PER_PAIR_INS = PER_PAIR_INS;
    sorter_mutator.PER_PAIR_DEL = PER_PAIR_DEL;
    sorter_mutator.PER_PAIR_SWAP = PER_PAIR_SWAP;

    SetMutFun([this](sorting_t & org, emp::Random & rnd) {
        // return sorter_mutator.Mutate(rnd, org);
        // std::cout << org.AsString() << org.GetSize() << std::endl;
        return sorter_mutator.Mutate(rnd, org);
        // std::cout << org.AsString() << org.GetSize() << std::endl << std::endl;
        // return 1.0;
    });

    for (size_t i = 0; i < START_POP_SIZE; ++i) {
        Inject(sorter_mutator.GenRandomBitSorter(*random_ptr));
    }    

}

template <>
void EcologyWorld<bit_t>::InitPop() {
    // Setup the mutation function.
    std::function<size_t(bit_t &, emp::Random &)> mut_fun =
        [this](bit_t & org, emp::Random & random) {
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

    // Build a random initial population
    for (uint32_t i = 0; i < START_POP_SIZE; i++) {
        bit_t next_org(N);
        for (uint32_t j = 0; j < N; j++) next_org[j] = random_ptr->P(0.5);
        Inject(next_org);
    }
}

template <>
void EcologyWorld<rv_t>::InitPop() {

    // Setup the mutation function.
    std::function<size_t(rv_t &, emp::Random &)> mut_fun =
        [this](rv_t & org, emp::Random & random) {
            for (uint32_t m = 0; m < N; m++) {
                org[m] += random_ptr->GetRandNormal(0, MUT_RATE); 
            }
            return 1;
        };
    SetMutFun( mut_fun );

    // Build a random initial population
    for (uint32_t i = 0; i < START_POP_SIZE; i++) {
        rv_t next_org(N);
        for (uint32_t j = 0; j < N; j++) next_org[j] = random_ptr->GetDouble(1);
        Inject(next_org);
    }
}

template <>
void EcologyWorld<gp_t>::InitPop() {

    // Setup the mutation function.
    std::function<size_t(gp_t &, emp::Random &)> mut_fun =
        [this](gp_t & org, emp::Random & r){
            int count = 0;
            for (size_t i = 0; i < org.GetSize(); ++i) {
                if (r.P(MUT_RATE)) {
                    org.RandomizeInst(i, r);
                    count++;
                }
                for (size_t j = 0; j < gp_t::base_t::INST_ARGS; j++) {
                    if (r.P(ARG_SUB_RATE)) {
                        org.genome.sequence[i].args[j] = r.GetUInt(org.CPU_SIZE);        
                        count++;
                    }
                }
                if (r.P(INS_MUT_RATE)) {
                    if (org.GetSize() < MAX_GENOME_SIZE) {
                        org.genome.sequence.insert(org.genome.sequence.begin() + (int)i, gp_t::gene_t());
                        org.RandomizeInst(i, r);
                        count++;
                    }
                }
                if (r.P(DEL_MUT_RATE)) {
                    if (org.GetSize() > 1) {
                        org.genome.sequence.erase(org.genome.sequence.begin() + (int)i);
                        count++;
                    }
                }

            }
            return count;
        }; 
    SetMutFun(mut_fun);


    emp::Random & random = GetRandom();
    for (size_t i = 0; i < START_POP_SIZE; i++) {
        gp_t cpu(&inst_set);
        cpu.PushRandom(random, GENOME_SIZE);
        Inject(cpu.GetGenome());
    }
}

template <>
void EcologyWorld<rv_t>::SetupRealValue() {
    // emp::vector<emp::vector<double>> grid_anchors(DIMENSIONS);
    // mid_point.resize(DIMENSIONS);
    // ubounds.resize(DIMENSIONS);
    // lbounds.resize(DIMENSIONS);
    // mut_std.resize(DIMENSIONS);

    // for (size_t dim = 0; dim < DIMENSIONS; ++dim) {
    //   grid_anchors[dim].resize(HINT_GRID_RES, 0.0);
    //   double l = eval_function->get_lbound(dim);
    //   double u = eval_function->get_ubound(dim);
    //   double step = (u - l)/(HINT_GRID_RES - 1);
    //   mid_point[dim] = l + ((u - l)/2);
    //   ubounds[dim] = u;
    //   lbounds[dim] = l;
    //   mut_std[dim] = (u - l)*MUTATION_STD;
    //   for (size_t i = 0; i < HINT_GRID_RES; ++i) {
    //     grid_anchors[dim][i] = l + (i*step);
    //   }
    // }
    // max_dist = CalcDist(ubounds, lbounds);

    // // Fill out key_points. (NOTE: currently not generic for more than 2 dimensions...)
    // score_ceil = eval_function->get_fitness_goptima();
    // score_floor = score_ceil;
    // for (size_t i = 0; i < HINT_GRID_RES; ++i) {
    //   for (size_t j = 0; j < HINT_GRID_RES; ++j) {
    //     key_points.emplace_back();
    //     key_points.back().resize(DIMENSIONS);
    //     key_points.back()[0] = grid_anchors[0][i];
    //     key_points.back()[1] = grid_anchors[1][j];
    //     const double key_point_score = eval_function->evaluate(key_points.back());
    //     score_floor = (key_point_score < score_floor) ? key_point_score : score_floor;
    //   }
    // }

    // for (size_t i = 0; i < key_points.size(); ++i) {
    //   phen.distances[i] = CalcDist(key_points[i], agent.genome); // Don't want to end up with divide-by-zero error.
    //   // phen.testcase_scores[i] = (phen.score < 0) ? (1+phen.distances[i]/max_dist)*phen.score : (1 - phen.distances[i]/max_dist)*phen.score;
    //   phen.testcase_scores[i] = (1 - (phen.distances[i]/max_dist)) * phen.transformed_score;
    // }

}

template <>
void EcologyWorld<bit_t>::SetupNK() {
    // Create landscape based on correct random seed
    landscape = emp::NKLandscape(N, K, *random_ptr);

    // Set-up fitness function (queries NK Landscape)
    evaluation_fun =
        [](bit_t & org){ return emp::Sum(org.phenotype); };

    if (SELECTION == (uint32_t) SELECTION_METHOD::SHARING) {
        SetSharedFitFun(evaluation_fun, [](bit_t & org1, bit_t & org2){return emp::calc_hamming_distance(org1.phenotype, org2.phenotype);}, SHARING_THRESHOLD, SHARING_ALPHA);
    } else {
        SetFitFun(evaluation_fun);
    }

    // Make a fitness function for each gene (set of sites that determine a fitness contribution)
    // in the bitstring

    for (size_t gene = 0; gene < N; gene++) {
        fit_set.push_back([this, gene](bit_t & org){return landscape.GetFitness(gene, org);});
    }
} 

template <>
double EcologyWorld<gp_t>::RunTestcase(gp_t & org, test_case_t testcase) {
    org.ResetHardware();
    for (size_t i = 0; i < testcase.first.size(); i++) {
        org.SetInput((int)i, testcase.first[i]);
    }
    org.Process(EVAL_TIME);
    double divisor = testcase.second;
    if (divisor == 0) {
        divisor = 1;
    }
    const std::unordered_map<int, double> & outputs = org.GetOutputs();
    int min_output = 100000;
    double result;
    for (auto out : outputs) {
        if (out.first < min_output) {
            min_output = out.first;
        }
    }

    if (outputs.size() != 0) {
        result = 1 / (std::abs(org.GetOutput(min_output) - testcase.second)/std::abs(divisor));
    } else {
        result = 0;
    }

    // emp_assert(std::abs(result) != INFINITY);
    if (result > 1000) {
        result = 1000;
    }
    return result;
}

template <>
void EcologyWorld<gp_t>::SetupProgramSynthesis() {
    testing_set.LoadTestcases(TESTSET_FILE);
    training_set.LoadTestcases(TRAINSET_FILE);

    inst_set = gp_t::inst_lib_t::DefaultInstLib();

    inst_set.AddInst("Null", [](gp_t::hardware_t & hw, const gp_t::inst_t & inst) {
        return;
    } , 0, "Null instruction");

    inst_set.AddInst("Dereference", [](gp_t::hardware_t & hw, const gp_t::inst_t & inst) {
        if (hw.regs[inst.args[0]] >= 0 && hw.regs[inst.args[0]] < hw.regs.size()) {
            hw.regs[inst.args[1]] = hw.regs[(size_t)hw.regs[inst.args[0]]];
        }
    } , 3, "Dereference register 0 and store result in register 1");

    for (size_t testcase = 0; testcase < N_TEST_CASES; ++testcase) {
        fit_set.push_back([testcase, this](gp_t & org) {
            return RunTestcase(org, training_set[testcase]);
        });
    }

    fun_calc_fitness_t goal_fun = [](gp_t & org) {
        return emp::Sum(org.phenotype);
    };

    if (SELECTION == (uint32_t) SELECTION_METHOD::SHARING) {
        SetSharedFitFun(goal_fun, [](gp_t & org1, gp_t & org2) {return (double)emp::calc_hamming_distance(org1.phenotype, org2.phenotype);}, SHARING_THRESHOLD, SHARING_ALPHA);
    } else {
        SetFitFun(goal_fun);
    }

    evaluation_fun = [this](gp_t & org) {
        double result = 0;
        for (size_t testcase = 0; testcase < N_TEST_CASES; ++testcase) {
            result += RunTestcase(org, testing_set[testcase]);
        }
        return result;
    };

}

template <>
void EcologyWorld<gp_t>::SetupLogic9() {
    inst_set = gp_t::inst_lib_t::DefaultInstLib();

    inst_set.AddInst("Null", [](gp_t::hardware_t & hw, const gp_t::inst_t & inst) {
        return;
    } , 0, "Null instruction");

    // Configure the tasks. 
    // Zero out task inputs.
    for (size_t i = 0; i < MAX_LOGIC_TASK_NUM_INPUTS; ++i) task_inputs[i] = 0;
    input_load_id = 0;

    // Add tasks to set.
    // NAND
    task_set.AddTask("NAND", [](taskset_t::Task & task, const std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS> & inputs) {
        const task_io_t a = inputs[0], b = inputs[1];
        task.solutions.emplace_back(~(a&b));
    }, "NAND task");
    // NOT
    task_set.AddTask("NOT", [](taskset_t::Task & task, const std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS> & inputs) {
        const task_io_t a = inputs[0], b = inputs[1];
        task.solutions.emplace_back(~a);
        task.solutions.emplace_back(~b);
    }, "NOT task");
    // ORN
    task_set.AddTask("ORN", [](taskset_t::Task & task, const std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS> & inputs) {
        const task_io_t a = inputs[0], b = inputs[1];
        task.solutions.emplace_back((a|(~b)));
        task.solutions.emplace_back((b|(~a)));
    }, "ORN task");
    // AND
    task_set.AddTask("AND", [](taskset_t::Task & task, const std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS> & inputs) {
        const task_io_t a = inputs[0], b = inputs[1];
        task.solutions.emplace_back(a&b);
    }, "AND task");
    // OR
    task_set.AddTask("OR", [](taskset_t::Task & task, const std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS> & inputs) {
        const task_io_t a = inputs[0], b = inputs[1];
        task.solutions.emplace_back(a|b);
    }, "OR task");
    // ANDN
    task_set.AddTask("ANDN", [](taskset_t::Task & task, const std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS> & inputs) {
        const task_io_t a = inputs[0], b = inputs[1];
        task.solutions.emplace_back((a&(~b)));
        task.solutions.emplace_back((b&(~a)));
    }, "ANDN task");
    // NOR
    task_set.AddTask("NOR", [](taskset_t::Task & task, const std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS> & inputs) {
        const task_io_t a = inputs[0], b = inputs[1];
        task.solutions.emplace_back(~(a|b));
    }, "NOR task");
    // XOR
    task_set.AddTask("XOR", [](taskset_t::Task & task, const std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS> & inputs) {
        const task_io_t a = inputs[0], b = inputs[1];
        task.solutions.emplace_back(a^b);
    }, "XOR task");
    // EQU
    task_set.AddTask("EQU", [](taskset_t::Task & task, const std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS> & inputs) {
        const task_io_t a = inputs[0], b = inputs[1];
        task.solutions.emplace_back(~(a^b));
    }, "EQU task");
    // ECHO
    task_set.AddTask("ECHO", [](taskset_t::Task & task, const std::array<task_io_t, MAX_LOGIC_TASK_NUM_INPUTS> & inputs) {
        const task_io_t a = inputs[0], b = inputs[1];
        task.solutions.emplace_back(a);
        task.solutions.emplace_back(b);
    }, "ECHO task");

    inst_set.AddInst("Nand", [](gp_t::hardware_t & hw, const gp_t::inst_t & inst) {
        hw.regs[inst.args[2]] = ~((int)hw.regs[inst.args[0]]&(int)hw.regs[inst.args[1]]);
    } , 3, "WM[ARG3]=~(WM[ARG1]&WM[ARG2])");

    // setup score
    evaluation_fun = [this](gp_t & org) {
        // Num unique tasks completed + (TOTAL TIME - COMPLETED TIME)
        ResetTasks();
        org.ResetHardware();

        for (size_t i = 0; i < MAX_LOGIC_TASK_NUM_INPUTS; i++) {
            org.SetInput((int)i, task_inputs[i]);
        }
        for (size_t t = 0; t < EVAL_TIME; ++t) {
            org.SingleProcess();

            int min_output = 100000;
            for (auto out : org.GetOutputs()) {
                if (out.first < min_output) {
                    min_output = out.first;
                }
            }

            if (min_output < 100000) {
                task_set.Submit((uint32_t)org.GetOutput(min_output), t);
                if (task_set.AllTasksCompleted()) {
                    break;
                }
            }
        }

        double score = task_set.GetUniqueTasksCompleted();
        if (task_set.GetAllTasksCompletedTime() > 0) score += (EVAL_TIME - task_set.GetAllTasksCompletedTime());
        return score;
    };
    
    // Add fitness functions to fit_set
    for (size_t task_num = 0; task_num < task_set.GetSize(); ++task_num) {
        fit_set.push_back([this, task_num](gp_t & org) {

            ResetTasks();
            org.ResetHardware();

            for (size_t i = 0; i < MAX_LOGIC_TASK_NUM_INPUTS; i++) {
                org.SetInput((int)i, task_inputs[i]);
            }
            for (size_t t = 0; t < EVAL_TIME; ++t) {
                org.SingleProcess();

                int min_output = 100000;
                for (auto out : org.GetOutputs()) {
                    if (out.first < min_output) {
                        min_output = out.first;
                    }
                }

                if (min_output < 100000) {
                    task_set.Submit((uint32_t)org.GetOutput(min_output), t);
                    if (task_set.GetTask(task_num).GetCompletionCnt() > 0) {
                        break;
                    }
                }
            }

            return (int) (task_set.GetTask(task_num).GetCompletionCnt() > 0);
        });
    }

    fit_set.push_back([this](gp_t & org) {
        ResetTasks();
        org.ResetHardware();        

        for (size_t i = 0; i < MAX_LOGIC_TASK_NUM_INPUTS; i++) {
            org.SetInput((int)i, task_inputs[i]);
        }

        for (size_t t = 0; t < EVAL_TIME; ++t) {
            org.SingleProcess();

            int min_output = 100000;
            for (auto out : org.GetOutputs()) {
                if (out.first < min_output) {
                    min_output = out.first;
                }
            }

            if (min_output < 100000) {
                task_set.Submit((uint32_t)org.GetOutput(min_output), t);
                if (task_set.AllTasksCompleted()) {
                    break;
                }
            }
        }

        double score = 0;
        if (task_set.GetAllTasksCompletedTime() > 0) score += (EVAL_TIME - task_set.GetAllTasksCompletedTime())/(double)EVAL_TIME;
        return score;
    });

    if (SELECTION == (uint32_t) SELECTION_METHOD::SHARING) {
        SetSharedFitFun(evaluation_fun, [](gp_t & org1, gp_t & org2) {return (double)emp::calc_hamming_distance(org1.phenotype, org2.phenotype);}, SHARING_THRESHOLD, SHARING_ALPHA);
    } else {
        SetFitFun(evaluation_fun);
    }

}

template <>
void EcologyWorld<sorting_t>::SetupSortingNetworks() {
    std::cout << "Max passes: " << (1 << NUM_BITS) << std::endl;
    evaluation_fun = [](sorting_t & org){
        return emp::Sum(org.phenotype);
    };

    if (SELECTION == (uint32_t) SELECTION_METHOD::SHARING) {
        SetSharedFitFun(evaluation_fun, [](sorting_t & org1, sorting_t & org2) {return (double)emp::calc_hamming_distance(org1.phenotype, org2.phenotype);}, SHARING_THRESHOLD, SHARING_ALPHA);
    } else {
        SetFitFun(evaluation_fun);
    }

    for (uint32_t i = 0; i < (uint32_t)(1 << NUM_BITS); i++) {
        fit_set.push_back([i](sorting_t & org){
            return org.TestSortable(i);
        });
    }
    
    fit_set.push_back([this](sorting_t & org){
        if (org.CountSortable(NUM_BITS) == (size_t)(1 << NUM_BITS)) {
            return ((double)MAX_NETWORK_SIZE - (double)org.GetSize())/(double)MAX_NETWORK_SIZE;
        }
        return 0.0;
    });

}