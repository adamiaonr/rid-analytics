#include <time.h>
#include <algorithm>
#include <string>

#include "argvparser.h"

#include "dataparser.h"
#include "rid-analytics.h"

#define MAX_ARRAY_SIZE      128

#define OPTION_SCN_FILE     (char *) "scn-file"
#define OPTION_DATA_DIR     (char *) "data-dir"
#define OPTION_VERBOSE      (char *) "verbose"
#define OPTION_TITLE        (char *) "test-title"
#define OPTION_CACHE_TIER   (char *) "cache-tier"

// output data in csv format for now
#define OUTPUT_FILE_EXT     (char *) "csv"

#define ALPHA_MAX           (double) 0.1
#define ALPHA_MIN           (double) 0.001
#define ALPHA_VALUES_SIZE   (int) 2
const double ALPHA_VALUES[] = {ALPHA_MIN, ALPHA_MAX};

const int _CONTENT_SOURCES[] = {1, 1, 1, 1};

#define PENALTY_TYPE_SIZE           2

#define FP_PROB             (char *) "fp_prob"
#define ALPHA               (char *) "alpha"
#define TIER_DEPTH          (char *) "tier_depth"
#define CONTENT_SOURCES     (char *) "content_sources"
#define DOMAINS             (char *) "domains"
#define LATENCIES           (char *) "latencies"

using namespace std;
using namespace CommandLineProcessing;

ArgvParser * create_argv_parser() {

    ArgvParser * cmds = new ArgvParser();

    cmds->setIntroductoryDescription("\n\nrid-analytics v0.1"\
        "\n\n *** RID-ANALYTICS SENSITIVITY ANALYSIS ***"\
        "\n\nrun simple analytical "\
        "evaluations on networks which use RIDs (basically Bloom Filters) for packet "\
        "forwarding. it computes probabilities of different outcomes and avg. latencies, "\
        "given a set of scenario characteristics (e.g. use of "\
        "caching, distribution of content sources, etc.).\ncreated by by adamiaonr@cmu.edu");

    cmds->setHelpOption("h", "help", "help page.");

    cmds->defineOption(
            OPTION_SCN_FILE,
            "path to .scn file which contains the input parameters for the "\
            "example (e.g. nr. of tiers, nr. of domains per tier, etc.)",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_DATA_DIR,
            "path to a directory to output .csv files w/ latency & outcome data.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_VERBOSE,
            "print stats during the model run. default: non-verbose.",
            ArgvParser::NoOptionAttribute);

    cmds->defineOptionAlternative(OPTION_VERBOSE, "v");

    cmds->defineOption(
            OPTION_TITLE,
            "common title to append to output file names (e.g. <TITLE>.outcomes.csv, <TITLE>.dot, etc.)",
            ArgvParser::OptionRequiresValue);

        cmds->defineOption(
            OPTION_CACHE_TIER,
            "first tier number at which content becomes visible.",
            ArgvParser::OptionRequiresValue);

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

    // parse the arguments with an ArgvParser object
    ArgvParser * cmds = create_argv_parser();

    // string for config file
    char scn_file[MAX_ARRAY_SIZE];
    // string for directory path where all files with output data will be
    char data_dir[MAX_ARRAY_SIZE];
    // string for output title (e.g. for .dot, .csv & other files)
    char title[MAX_ARRAY_SIZE];
    // cache location define as input argument
    char cache_tier_str[MAX_ARRAY_SIZE];
    int cache_tier = 0;
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

            strncpy(scn_file, (char *) cmds->optionValue(OPTION_SCN_FILE).c_str(), MAX_ARRAY_SIZE);

        } else {

            fprintf(stderr, "no .scn file path specified. use "\
                "option -h for help.\n");

            delete cmds;
            return -1;
        }

        if (cmds->foundOption(OPTION_DATA_DIR)) {

            strncpy(data_dir, (char *) cmds->optionValue(OPTION_DATA_DIR).c_str(), MAX_ARRAY_SIZE);
        
        } else {

            fprintf(stderr, "no data directory path specified. use "\
                "option -h for help.\n");

            delete cmds;
            return -1;
        }

        if (cmds->foundOption(OPTION_TITLE)) {

            strncpy(title, (char *) cmds->optionValue(OPTION_TITLE).c_str(), MAX_ARRAY_SIZE);
        
        } else {

            fprintf(stderr, "no title specified for output files. use "\
                "option -h for help.\n");

            delete cmds;
            return -1;
        }

        if (cmds->foundOption(OPTION_CACHE_TIER)) {

            strncpy(cache_tier_str, (char *) cmds->optionValue(OPTION_CACHE_TIER).c_str(), MAX_ARRAY_SIZE);
            cache_tier = atoi(cache_tier_str);
        
        } else {

            fprintf(stderr, "no title specified for output files. use "\
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

    // extract the sensitivity analysis from the .scn file
    DataParser * scn_parser = new DataParser(scn_file);

    // nr. of tiers
    int tier_depth = 0;
    scn_parser->get_int_property_value(TIER_DEPTH, tier_depth);

    // array of fp probabilities to test
    double * fp_prob = (double *) calloc(MAX_ARRAY_SIZE, sizeof(double));
    scn_parser->get_double_property_array(FP_PROB, fp_prob);
    int fp_prob_size = get_array_size(fp_prob);

    // array of alphas to test
    double * alpha = (double *) calloc(MAX_ARRAY_SIZE, sizeof(double));
    scn_parser->get_double_property_array(ALPHA, alpha);
    int alpha_size = get_array_size(alpha);

    // array of latencies per tier
    double * latencies = (double *) calloc(tier_depth, sizeof(double));
    scn_parser->get_double_property_array(LATENCIES, latencies);

    // array of content sources per tier (this is calculated from the 
    // cache_tier input arg)
    int * content_sources = (int *) calloc(tier_depth, sizeof(int));
    for (int i = (tier_depth - 1); (i + 1) > cache_tier; i--)
        content_sources[i] = _CONTENT_SOURCES[i];

    // ... and domains per tier
    int * domains = (int *) calloc(tier_depth, sizeof(int));
    scn_parser->get_int_property_array(DOMAINS, domains);

    // let the games begin...

    // arrays which will keep the test values during the analysis
    double * _fp_prob = (double *) calloc(tier_depth, sizeof(double));
    double * _alpha = (double *) calloc(tier_depth, sizeof(double));

    // arrays of FILE * are the best interface for the job, even if only 1 FILE * 
    // is used. btw, the interfaces are defined in rid-analytics.h.
    FILE ** latencies_files = (FILE **) calloc(1, sizeof(FILE *));
    FILE ** outcomes_files = (FILE **) calloc(1, sizeof(FILE *));

    // open the output .csv files in "a"ppend mode
    latencies_files[0] = fopen(
        std::string(std::string(data_dir) + "/latencies." + OUTPUT_FILE_EXT).c_str(), 
        "a");

    outcomes_files[0] = fopen(
        std::string(std::string(data_dir) + "/outcomes." + OUTPUT_FILE_EXT).c_str(), 
        "a");

    double avg_latency = 0.0;
    clock_t begin, end;
    unsigned long int nr_tests = 0;

    // struct for input parameters
    RIDAnalytics::rid_analytics_inputs input_params;

    // these input parameters won't change during the run
    input_params.tier_depth = tier_depth;
    input_params.content_sources = content_sources;
    input_params.domains = domains;
    input_params.latencies = latencies;
    input_params.title = title;

    // keep track of execution time
    begin = clock();

    // the idea now is to have all tiers with the same value (both _fp_prob and 
    // alpha)
    for (int i = 0; i < fp_prob_size; i++) {

        for (int j = 0; j < tier_depth; j++) {

            _fp_prob[j] = fp_prob[i];
        }

        // update the fp_prob input param
        input_params.fp_prob = _fp_prob;

        for (int k = 0; k < alpha_size; k++) {

            for (int m = 0; m < tier_depth; m++) {

                _alpha[m] = alpha[k];                
            }

            // update the alpha input param
            input_params.alpha = _alpha;

            // for the value returned by run_model()
            avg_latency = 0.0;

            // RUN #1 : FEEDBACK penalty type
            input_params.fp_resolution_tech = PENALTY_TYPE_FEEDBACK;

            RIDAnalytics::run_model(
                avg_latency,
                input_params,
                (verbose | MODE_SAVEOUTCOMES),
                outcomes_files,
                1,
                title,
                std::string(data_dir));

            // write to the .csv files, for each _fp_prob[tier]
            // value
            // format should be:
            // [fp@<_fp_tier>][_fp_prob[<_fp_tier>]][a@<_alpha_tier>][_alpha[_alpha_tier]][title][penalty_type][avg_latency]
            for (int t = 0; t < tier_depth; t++) {

                fprintf(
                    latencies_files[0], 
                    "%-.8E,%-.8E,feedback,%-.8E\n", _fp_prob[t], _alpha[t], avg_latency);
            }

            // RUN #2 : FALLBACK penalty type
            avg_latency = 0.0;
            input_params.fp_resolution_tech = PENALTY_TYPE_FALLBACK;

            RIDAnalytics::run_model(
                avg_latency,
                input_params,
                (verbose | MODE_SAVEGRAPH | MODE_SAVEOUTCOMES),
                outcomes_files,
                1,
                title,
                std::string(data_dir));

            for (int t = 0; t < tier_depth; t++) {

                fprintf(
                    latencies_files[0], 
                    "%-.8E,%-.8E,fallback,%-.8E\n", _fp_prob[t], _alpha[t], avg_latency);
            }

            nr_tests++;
        }
    }

    // close files
    fclose(latencies_files[0]);
    fclose(outcomes_files[0]);

    end = clock();

    printf("[EXECUTION TIME : %-.8f sec]\n", (double)(end - begin) / CLOCKS_PER_SEC);
    printf("[# OF RUNS : %-8ld]\n", nr_tests);

    // clean everything up...
    delete scn_parser;
    delete cmds;
    
    free(fp_prob);
    free(alpha);
    free(_fp_prob);
    free(_alpha);
    free(latencies);
    free(content_sources);
    free(domains);

    free(latencies_files);
    free(outcomes_files);

    return 0;
}