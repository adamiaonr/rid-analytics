#include <time.h>

#include "argvparser.h"

#include "dataparser.h"
#include "rid-analytics.h"

#define MAX_ARRAY_SIZE      128

#define OPTION_SCN_FILE     (char *) "scn-file"
#define OPTION_DATA_DIR     (char *) "data-dir"
#define OPTION_VERBOSE      (char *) "verbose"

#define FP_PROB             (char *) "fp_prob"
#define O_OPTIMISTIC        (char *) "o_optimistic"
#define ORIGIN_LEVEL        (char *) "origin_level"
#define LEVEL_DEPTH         (char *) "level_depth"

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
    bool verbose = false;

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
            verbose = true;
        }
    }

    if (result == ArgvParser::ParserHelpRequested) {

        delete cmds;
        return -1;
    }

    // 2) extract the sensitivity analysis from the .scn file
    DataParser * scn_parser = new DataParser(scn_file);

    // nr. of levels
    int level_depth = 0;
    scn_parser->get_int_property_value(LEVEL_DEPTH, level_depth);

    // array of fp probability values (each value will be tested for all levels) 
    double * fp_prob = (double *) calloc(MAX_ARRAY_SIZE, sizeof(double));
    scn_parser->get_double_property_array(FP_PROB, fp_prob);
    int fp_prob_size = get_array_size(fp_prob);

    // array of o-optimistic values
    double * o_optimistic = (double *) calloc(MAX_ARRAY_SIZE, sizeof(double));
    scn_parser->get_double_property_array(O_OPTIMISTIC, o_optimistic);
    int o_optimistic_size = get_array_size(o_optimistic);

    // 'levels' at which origin and cache are located
    int origin_level = 0;
    scn_parser->get_int_property_value(ORIGIN_LEVEL, origin_level);
    int cache_level = 3;

    // // fp resolution mechanisms
    // int fp_resolution[2] = {PENALTY_TYPE_FEEDBACK, PENALTY_TYPE_FALLBACK};
    // int fp_resolution_size = sizeof(fp_resolution) / sizeof(int);

    // let the games begin...

    // arrays which will keep the test values during the analysis
    double * _fp_prob = (double *) calloc(level_depth, sizeof(double));
    double * _o_optimistic = (double *) calloc(level_depth, sizeof(double));

    // arrays of FILE * (this can only end well... LOL)
    FILE ** fp_prob_file = (FILE **) calloc(level_depth, sizeof(FILE *));
    FILE ** fp_prob_outcomes_file = (FILE **) calloc(level_depth, sizeof(FILE *));

    FILE ** o_optimistic_file = (FILE **) calloc(level_depth, sizeof(FILE *));
    FILE ** o_optimistic_outcomes_file = (FILE **) calloc(level_depth, sizeof(FILE *));

    // open the .csv files in "a"ppend mode
    for (int i = 0; i < level_depth; i++) {

        fp_prob_file[i] = fopen(
            std::string(std::string(data_dir) + "/fp." + std::to_string(i + 1) + ".csv").c_str(), 
            "a");

        fp_prob_outcomes_file[i] = fopen(
            std::string(std::string(data_dir) + "/fp." + std::to_string(i + 1) + ".outcomes.csv").c_str(), 
            "a");

        o_optimistic_file[i] = fopen(
            std::string(std::string(data_dir) + "/op." + std::to_string(i + 1) + ".csv").c_str(), 
            "a");

        o_optimistic_outcomes_file[i] = fopen(
            std::string(std::string(data_dir) + "/op." + std::to_string(i + 1) + ".outcomes.csv").c_str(), 
            "a");
    }

    double avg_latency = 0.0;
    clock_t begin, end;
    unsigned long int nr_tests = 0;

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

                    for (int a = 0; a < o_optimistic_size; a++) {
                        _o_optimistic[0] = o_optimistic[a];

                        for (int b = 0; b < o_optimistic_size; b++) {
                            _o_optimistic[1] = o_optimistic[b];

                            for (int c = 0; c < o_optimistic_size; c++) {
                                _o_optimistic[2] = o_optimistic[c];

                                for (int d = 0; d < o_optimistic_size; d++) {
                                    _o_optimistic[3] = o_optimistic[d];

                                    avg_latency = 0.0;

                                    RIDAnalytics::run_model(
                                        avg_latency,
                                        _fp_prob, 
                                        _o_optimistic, 
                                        level_depth,
                                        cache_level, 
                                        origin_level,
                                        PENALTY_TYPE_FEEDBACK,
                                        verbose,
                                        false,
                                        false,
                                        false,
                                        NULL,
                                        NULL,
                                        std::string(data_dir));

                                    // write to the .csv files
                                    for (int f = 0; f < level_depth; f++) {

                                        fprintf(
                                            fp_prob_file[f], 
                                            "%-.8E,%-.8E,feedback\n", _fp_prob[f], avg_latency);
                                        fprintf(
                                            o_optimistic_file[f], 
                                            "%-.8E,%-.8E,feedback\n", _o_optimistic[f], avg_latency);
                                    }

                                    avg_latency = 0.0;

                                    RIDAnalytics::run_model(
                                        avg_latency,
                                        _fp_prob, 
                                        _o_optimistic, 
                                        level_depth,
                                        cache_level, 
                                        origin_level,
                                        PENALTY_TYPE_FALLBACK,
                                        verbose,
                                        false,
                                        false,
                                        true,
                                        fp_prob_outcomes_file,
                                        o_optimistic_outcomes_file,
                                        std::string(data_dir));

                                    // write to the .csv files
                                    for (int f = 0; f < level_depth; f++) {

                                        fprintf(
                                            fp_prob_file[f], 
                                            "%-.8E,%-.8E,fallback\n", _fp_prob[f], avg_latency);
                                        fprintf(
                                            o_optimistic_file[f], 
                                            "%-.8E,%-.8E,fallback\n", _o_optimistic[f], avg_latency);
                                    
                                        // fflush(fp_prob_file[f]);
                                        // fflush(o_optimistic_file[f]);
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
    for (int f = 0; f < level_depth; f++) {
        fclose(fp_prob_file[f]);
        fclose(fp_prob_outcomes_file[f]);

        fclose(o_optimistic_file[f]);
        fclose(o_optimistic_outcomes_file[f]);
    }

    end = clock();

    printf("[EXECUTION TIME : %-.8f sec]\n", (double)(end - begin) / CLOCKS_PER_SEC);
    printf("[# OF RUNS : %-8ld]\n", nr_tests);

    // clean everything up...
    delete scn_parser;
    delete cmds;
    free(_fp_prob);
    free(_o_optimistic);
    free(fp_prob_file);
    free(fp_prob_outcomes_file);
    free(o_optimistic_file);
    free(o_optimistic_outcomes_file);

    return 0;
}