#include <time.h>

#include "argvparser.h"

#include "dataparser.h"
#include "rid-analytics.h"

#define MAX_ARRAY_SIZE      128

#define OPTION_SCN_FILE     (char *) "scn-file"
#define OPTION_DATA_DIR     (char *) "data-dir"
#define OPTION_VERBOSE      (char *) "verbose"

#define FP_PROB             (char *) "fp_prob"
#define ALPHA               (char *) "alpha"
#define ORIGIN_TIER         (char *) "origin_tier"
#define CACHE_TIER          (char *) "cache_tier"
#define TIER_DEPTH          (char *) "tier_depth"

using namespace std;
using namespace CommandLineProcessing;

ArgvParser * create_argv_parser() {

    ArgvParser * cmds = new ArgvParser();

    cmds->setIntroductoryDescription("\n\nrid-analytics v0.1\n\ntool to run simple analytical "\
        "evaluations on networks which use RIDs (i.e. Bloom Filters) for packet "\
        "forwarding. it computes probability distributions for latencies and "\
        "avg. latencies/ relaying penalties given a set of scenario characteristics (e.g. use of "\
        "caching, distribution of content sources, etc.).\nby adamiaonr@cmu.edu");

    cmds->setHelpOption("h", "help", "help page.");

    cmds->defineOption(
            OPTION_SCN_FILE,
            "path to .scn file which contains the input info for the sensitivity analysis",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_DATA_DIR,
            "path to a directory to output .csv files w/ latency data.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_VERBOSE,
            "print stats during the model run. default: not verbose.",
            ArgvParser::NoOptionAttribute);
    cmds->defineOptionAlternative(OPTION_VERBOSE, "v");

    return cmds;
}

int get_array_size(double * array) {

    int size = 0;

    for (int i = 0; i < MAX_ARRAY_SIZE; i++) {

        if (!(array[i] > 0.0))
            break;

        size++;
    }

    return size;
}

int main (int argc, char **argv) {

    // 1) first order of business: parse the arguments with an ArgvParser
    ArgvParser * cmds = create_argv_parser();

    // string for config file
    char * scn_file;
    // string for directory path where all files with output data will be
    char * data_dir;
    // keeps track of verbose mode
    unsigned int verbose = 0;

    // parse() takes the arguments to main() and parses them according to 
    // ArgvParser rules
    int result = cmds->parse(argc, argv);

    // if something went wrong: show help option
    if (result != ArgvParser::NoParserError) {

        fprintf(stderr, "%s\n", cmds->parseErrorDescription(result).c_str());

        if (result != ArgvParser::ParserHelpRequested) {
            fprintf(stderr, "use option -h for help.\n");
        }

        delete cmds;
        return -1;

    } else {

        // otherwise, check for the different OPTION_
        if (cmds->foundOption(OPTION_SCN_FILE)) {

            scn_file = (char *) cmds->optionValue(OPTION_SCN_FILE).c_str();

        } else {

            fprintf(stderr, "no .scn file path specified. use "\
                "option -h for help.\n");

            delete cmds;
            return -1;
        }

        if (cmds->foundOption(OPTION_DATA_DIR)) {

            data_dir = (char *) cmds->optionValue(OPTION_DATA_DIR).c_str();
        
        } else {

            fprintf(stderr, "no data directory path specified. use "\
                "option -h for help.\n");

            delete cmds;
            return -1;
        }

        if (cmds->foundOption(OPTION_VERBOSE)) {
            verbose = MODE_VERBOSE;
        }
    }

    if (result == ArgvParser::ParserHelpRequested) {

        delete cmds;
        return -1;
    }

    // 2) extract the sensitivity analysis from the .scn file
    DataParser * scn_parser = new DataParser(scn_file);

    // nr. of tiers
    int tier_depth = 0;
    scn_parser->get_int_property_value(TIER_DEPTH, tier_depth);

    // array of fp probability values (each value will be tested for all tiers) 
    double * fp_prob = (double *) calloc(MAX_ARRAY_SIZE, sizeof(double));
    scn_parser->get_double_property_array(FP_PROB, fp_prob);
    int fp_prob_size = get_array_size(fp_prob);

    // array of alpha values
    double * alpha = (double *) calloc(MAX_ARRAY_SIZE, sizeof(double));
    scn_parser->get_double_property_array(ALPHA, alpha);
    int alpha_size = get_array_size(alpha);

    // 'tiers' at which origin and cache are located
    int origin_tier = 3;
    scn_parser->get_int_property_value(ORIGIN_TIER, origin_tier);
    int cache_tier = 2;
    scn_parser->get_int_property_value(CACHE_TIER, cache_tier);

    // // fp resolution mechanisms
    // int fp_resolution[2] = {PENALTY_TYPE_FEEDBACK, PENALTY_TYPE_FALLBACK};
    // int fp_resolution_size = sizeof(fp_resolution) / sizeof(int);

    // let the games begin...

    // arrays which will keep the test values during the analysis
    double * _fp_prob = (double *) calloc(tier_depth, sizeof(double));
    double * _alpha = (double *) calloc(tier_depth, sizeof(double));

    // arrays of FILE * (this can only end well... LOL)
    FILE ** fp_prob_file = (FILE **) calloc(tier_depth, sizeof(FILE *));
    FILE ** fp_prob_outcomes_file = (FILE **) calloc(tier_depth, sizeof(FILE *));
    // FILE ** alpha_file = (FILE **) calloc(tier_depth, sizeof(FILE *));
    // FILE ** alpha_outcomes_file = (FILE **) calloc(tier_depth, sizeof(FILE *));

    // open the .csv files in "a"ppend mode
    for (int i = 0; i < tier_depth; i++) {

        fp_prob_file[i] = fopen(
            std::string(std::string(data_dir) + "/fp." + std::to_string(i + 1) + ".csv").c_str(), 
            "a");

        fp_prob_outcomes_file[i] = fopen(
            std::string(std::string(data_dir) + "/fp." + std::to_string(i + 1) + ".outcomes.csv").c_str(), 
            "a");

        // alpha_file[i] = fopen(
        //     std::string(std::string(data_dir) + "/op." + std::to_string(i + 1) + ".csv").c_str(), 
        //     "a");

        // alpha_outcomes_file[i] = fopen(
        //     std::string(std::string(data_dir) + "/op." + std::to_string(i + 1) + ".outcomes.csv").c_str(), 
        //     "a");
    }

    double avg_latency = 0.0;
    clock_t begin, end;
    unsigned long int nr_tests = 0;

    // struct for input parameters
    RIDAnalytics::rid_analytics_inputs input_params;
    // these input parameters won't change during the run
    input_params.tier_depth = tier_depth;
    input_params.cache_tier = cache_tier;
    input_params.origin_tier = origin_tier;

    // keep track of execution time
    begin = clock();

    for (int i = 0; i < fp_prob_size; i++) {
        _fp_prob[0] = fp_prob[i];

        for (int j = 0; j < fp_prob_size; j++) {
            _fp_prob[1] = fp_prob[j];

            for (int k = 0; k < fp_prob_size; k++) {
                _fp_prob[2] = fp_prob[k];

                for (int l = 0; l < fp_prob_size; l++) {
                    _fp_prob[3] = fp_prob[l];

                    for (int a = 0; a < alpha_size; a++) {
                        _alpha[0] = alpha[a];

                        for (int b = 0; b < alpha_size; b++) {
                            _alpha[1] = alpha[b];

                            for (int c = 0; c < alpha_size; c++) {
                                _alpha[2] = alpha[c];

                                for (int d = 0; d < alpha_size; d++) {
                                    _alpha[3] = alpha[d];

                                    avg_latency = 0.0;

                                    // update the model inputs
                                    input_params.fp_prob = _fp_prob;
                                    input_params.alpha = _alpha;
                                    input_params.fp_resolution_tech = PENALTY_TYPE_FEEDBACK;

                                    RIDAnalytics::run_model(
                                        avg_latency,
                                        input_params,
                                        verbose,
                                        NULL,
                                        std::string(data_dir));

                                    // write to the .csv files
                                    for (int f = 0; f < tier_depth; f++) {

                                        fprintf(
                                            fp_prob_file[f], 
                                            "%-.8E,%-.8E,feedback\n", _fp_prob[f], avg_latency);
                                        // fprintf(
                                        //     alpha_file[f], 
                                        //     "%-.8E,%-.8E,feedback\n", _alpha[f], avg_latency);
                                    }

                                    avg_latency = 0.0;
                                    // IMPORTANT: change penalty type to 
                                    // fallback
                                    input_params.fp_resolution_tech = PENALTY_TYPE_FALLBACK;

                                    RIDAnalytics::run_model(
                                        avg_latency,
                                        input_params,
                                        (verbose | MODE_SAVEOUTCOMES),
                                        fp_prob_outcomes_file,
                                        std::string(data_dir));

                                    // write to the .csv files
                                    for (int f = 0; f < tier_depth; f++) {

                                        fprintf(
                                            fp_prob_file[f], 
                                            "%-.8E,%-.8E,fallback\n", _fp_prob[f], avg_latency);
                                        // fprintf(
                                        //     alpha_file[f], 
                                        //     "%-.8E,%-.8E,fallback\n", _alpha[f], avg_latency);
                                    
                                        // fflush(fp_prob_file[f]);
                                        // fflush(alpha_file[f]);
                                    }

                                    nr_tests++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // close files
    for (int f = 0; f < tier_depth; f++) {
        fclose(fp_prob_file[f]);
        fclose(fp_prob_outcomes_file[f]);

        // fclose(alpha_file[f]);
        // fclose(alpha_outcomes_file[f]);
    }

    end = clock();

    printf("[EXECUTION TIME : %-.8f sec]\n", (double)(end - begin) / CLOCKS_PER_SEC);
    printf("[# OF RUNS : %-8ld]\n", nr_tests);

    // clean everything up...
    delete scn_parser;
    delete cmds;
    free(_fp_prob);
    free(_alpha);
    free(fp_prob_file);
    free(fp_prob_outcomes_file);
    // free(alpha_file);
    // free(alpha_outcomes_file);

    return 0;
}