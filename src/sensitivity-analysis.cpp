#include "argvparser.h"

#include "dataparser.h"
#include "rid-analytics.h"

#define MAX_ARRAY_SIZE      128

#define OPTION_SCN_FILE     (char *) "scn-file"
#define OPTION_DATA_DIR     (char *) "data-dir"

#define FP_PROB             (char *) "fp_prob"
#define O_OPTIMISTIC        (char *) "o_optimistic"
#define ORIGIN_LEVEL        (char *) "origin_level"
#define LEVEL_DEPTH         (char *) "level_depth"

#define PROB_SIZE           5
#define O_OPTIMISTIC_SIZE   4

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

    return cmds;
}

int main (int argc, char **argv) {

    ArgvParser * cmds = create_argv_parser();

    char * scn_file;
    char * data_dir;
    // trigger argument parsing with parse()
    int result = cmds->parse(argc, argv);

    if (result != ArgvParser::NoParserError) {

        // something went wrong. show help option.
        fprintf(stderr, "%s\n", cmds->parseErrorDescription(result).c_str());

        if (result != ArgvParser::ParserHelpRequested) {
            fprintf(stderr, "use option -h for help.\n");
        }

        return -1;

    } else {

        if (cmds->foundOption(OPTION_SCN_FILE)) {

            scn_file = (char *) cmds->optionValue(OPTION_SCN_FILE).c_str();

        } else {

            fprintf(stderr, "no .scn file path specified. use "\
                "option -h for help.\n");

            return -1;
        }

        if (cmds->foundOption(OPTION_DATA_DIR)) {

            data_dir = (char *) cmds->optionValue(OPTION_DATA_DIR).c_str();
        
        } else {

            fprintf(stderr, "no data directory path specified. use "\
                "option -h for help.\n");

            return -1;
        }
    }

    if (result == ArgvParser::ParserHelpRequested) {
        return -1;
    }

    // 1) extract the sensitivity analysis from the .scn file
    DataParser * scn_parser = new DataParser(scn_file);

    // 1.1) nr. of levels
    int level_depth = 0;
    scn_parser->get_int_property_value(LEVEL_DEPTH, level_depth);

    // 1.2) fp probabilities
    double fp_prob[MAX_ARRAY_SIZE];
    scn_parser->get_double_property_array(FP_PROB, fp_prob);

    // 1.3) o-optimistic values
    double o_optimistic[MAX_ARRAY_SIZE];
    scn_parser->get_double_property_array(O_OPTIMISTIC, o_optimistic);

    // 1.4) origin and cache levels
    int origin_level = 0;
    scn_parser->get_int_property_value(ORIGIN_LEVEL, origin_level);
    int cache_level = 2;

    // 1.5) fp resolution mechanisms
    //int fp_resolution_tech[2] = {PENALTY_TYPE_FEEDBACK, PENALTY_TYPE_FALLBACK};

    // 2) let the games begin...
    char float_str[16];
    double * _fp_prob = (double *) calloc(level_depth, sizeof(double));
    double * _o_optimistic = (double *) calloc(level_depth, sizeof(double));

    // 2.1) arrays of filestreams... LOL
    ofstream fp_prob_filestreams[4];
    ofstream o_optimistic_filestreams[4];

    // for (int i = 0; i < level_depth; i++) {
    //     fp_prob_filestreams[i].open(
    //         std::string(std::string(data_dir) + "/fp." + std::to_string(i + 1) + ".csv").c_str(), 
    //         std::fstream::in | std::fstream::out | std::fstream::app);
    // }

    // for (int i = 0; i < level_depth; i++) {
    //     o_optimistic_filestreams[i].open(
    //         std::string(std::string(data_dir) + "/op." + std::to_string(i + 1) + ".csv").c_str(), 
    //         std::fstream::in | std::fstream::out | std::fstream::app);
    // }

    double avg_latency = 0.0;

    for (int i = 0; i < PROB_SIZE; i++) {
        _fp_prob[0] = fp_prob[i];

        for (int j = 0; j < PROB_SIZE; j++) {
            _fp_prob[1] = fp_prob[j];

            for (int k = 0; k < PROB_SIZE; k++) {
                _fp_prob[2] = fp_prob[k];

                for (int l = 0; l < PROB_SIZE; l++) {
                    _fp_prob[3] = fp_prob[l];

                    for (int a = 0; a < O_OPTIMISTIC_SIZE; a++) {
                        _o_optimistic[0] = o_optimistic[a];

                        for (int b = 0; b < O_OPTIMISTIC_SIZE; b++) {
                            _o_optimistic[1] = o_optimistic[b];

                            for (int c = 0; c < O_OPTIMISTIC_SIZE; c++) {
                                _o_optimistic[2] = o_optimistic[c];

                                for (int d = 0; d < O_OPTIMISTIC_SIZE; d++) {
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
                                        false,
                                        false,
                                        std::string(data_dir));

                                    // for (int f = 0; f < level_depth; f++) {
                                        
                                    //     memset(float_str, 0, 16 * sizeof(char));
                                    //     snprintf(float_str, 16, "%-.8E", _fp_prob[f]);
                                    //     fp_prob_filestreams[f] << std::string(float_str) + ",";

                                    //     memset(float_str, 0, 16 * sizeof(char));
                                    //     snprintf(float_str, 16, "%-.8E", avg_latency);
                                    //     fp_prob_filestreams[f] << std::string(float_str) + ",feedback\n";
                                    // }

                                    // for (int f = 0; f < level_depth; f++) {
                                        
                                    //     memset(float_str, 0, 16 * sizeof(char));
                                    //     snprintf(float_str, 16, "%-.8E", _o_optimistic[f]);
                                    //     o_optimistic_filestreams[f] << std::string(float_str) + ",";

                                    //     memset(float_str, 0, 16 * sizeof(char));
                                    //     snprintf(float_str, 16, "%-.8E", avg_latency);
                                    //     o_optimistic_filestreams[f] << std::string(float_str) + ",feedback\n";
                                    // }

                                    avg_latency = 0.0;

                                    RIDAnalytics::run_model(
                                        avg_latency,
                                        _fp_prob, 
                                        _o_optimistic, 
                                        level_depth,
                                        cache_level, 
                                        origin_level,
                                        PENALTY_TYPE_FALLBACK,
                                        false,
                                        false,
                                        std::string(data_dir));

                                    // for (int f = 0; f < level_depth; f++) {
                                            
                                    //     memset(float_str, 0, 16 * sizeof(char));
                                    //     snprintf(float_str, 16, "%-.8E", _fp_prob[f]);
                                    //     fp_prob_filestreams[f] << std::string(float_str) + ",";

                                    //     memset(float_str, 0, 16 * sizeof(char));
                                    //     snprintf(float_str, 16, "%-.8E", avg_latency);
                                    //     fp_prob_filestreams[f] << std::string(float_str) + ",fallback\n";
                                    // }

                                    // for (int f = 0; f < level_depth; f++) {
                                        
                                    //     memset(float_str, 0, 16 * sizeof(char));
                                    //     snprintf(float_str, 16, "%-.8E", _o_optimistic[f]);
                                    //     o_optimistic_filestreams[f] << std::string(float_str) + ",";

                                    //     memset(float_str, 0, 16 * sizeof(char));
                                    //     snprintf(float_str, 16, "%-.8E", avg_latency);
                                    //     o_optimistic_filestreams[f] << std::string(float_str) + ",fallback\n";
                                    // }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}