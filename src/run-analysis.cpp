#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "argvparser.h"
#include "rid-analytics.h"

#define MAX_ARRAY_SIZE      256

#define OPTION_SCN_FILE             (char *) "scn-file"
#define OPTION_DATA_DIR             (char *) "data-dir"
#define OPTION_OUTPUT_LABEL         (char *) "output-label"
#define OPTION_VERBOSE              (char *) "verbose"

// eval parameters can be specified directly through the CLI and overwrite any 
// parameters set in .scn files
#define OPTION_REQUEST_SIZE         (char *) "request-size"
#define OPTION_BF_SIZE              (char *) "bf-size"
#define OPTION_MM_MODE              (char *) "mm-mode"
#define OPTION_EH_MODE              (char *) "eh-mode"
#define OPTION_RESOLV_MODE          (char *) "resolution-mode"
#define OPTION_ORIGIN_SERVER        (char *) "origin-server"
#define OPTION_START_ROUTER         (char *) "start-router"

using namespace std;
using namespace CommandLineProcessing;

ArgvParser * create_argv_parser() {

    ArgvParser * cmds = new ArgvParser();

    cmds->setIntroductoryDescription("\n\nrid-analytics v2.0\n\na tool to \
analyze the behavior of Bloom-filter based forwarding networks. it computes the \
probability of different delivery states for all possible request paths in a \
network topology passed as input.\nby adamiaonr@cmu.edu");

    cmds->setHelpOption("h", "help", "help page.");

    cmds->defineOption(
            OPTION_SCN_FILE,
            "path to .scn file with scenario info for analysis",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_DATA_DIR,
            "path to directory on which to save .csv files w/ results",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_OUTPUT_LABEL,
            "label to add to output .tsv files (will be created in <data-dir>)",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_BF_SIZE,
            "size of Bloom filter (in bits)",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_REQUEST_SIZE,
            "size of request, in # of URL components",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_MM_MODE,
            "multiple match resolution mode. 0 for 'FLOOD', 1 for 'RANDOM', 2 for 'FALLBACK'. default is 'FLOOD'.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_EH_MODE,
            "incorrect delivery handling mode. 0 for 'FEEDBACK', 1 for 'FALLBACK'. default is 'FEEDBACK'.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_RESOLV_MODE,
            "enable/disable error resolution. 0 for 'DISABLE', 1 for 'ENABLE'. default is 'DISABLE'.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_ORIGIN_SERVER,
            "id of origin server. default is '1'.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_START_ROUTER,
            "id of starting router. default is '0'.",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_VERBOSE,
            "print stats during the model run. default: not verbose",
            ArgvParser::NoOptionAttribute);
    cmds->defineOptionAlternative(OPTION_VERBOSE, "v");

    return cmds;
}

int main (int argc, char **argv) {

    // parse the arguments with an ArgvParser
    ArgvParser * cmds = create_argv_parser();

    // .scn config file path
    char scn_file[MAX_ARRAY_SIZE];
    // dir where all .csv files with output data will be saved
    char data_dir[MAX_ARRAY_SIZE];
    // output file name for this run
    char output_label[MAX_ARRAY_SIZE];
    // verbose mode (no verbosity by default)
    bool verbose = false;
    // evaluation parameters
    // request size (also largest possible |F|) in # of URL parameters
    int request_size = 0;       
    // bloom filter size (in bit)
    int bf_size = 0;
    // multiple match resolve mode
    int mm_mode = 0;
    // incorrect delivery handling mode
    int eh_mode = 0;
    int resolv_mode = 0;
    // origin server location. default is '1'
    char origin_server[MAX_ARRAY_SIZE];
    strncpy(origin_server, "1", MAX_ARRAY_SIZE);
    // start router id. default is '0'
    char start_router[MAX_ARRAY_SIZE];
    strncpy(start_router, "0", MAX_ARRAY_SIZE);

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

        if (cmds->foundOption(OPTION_OUTPUT_LABEL)) {

            strncpy(output_label, (char *) cmds->optionValue(OPTION_OUTPUT_LABEL).c_str(), MAX_ARRAY_SIZE);

        } else {

            fprintf(stderr, "no output .csv file name specified. use "\
                "option -h for help.\n");

            delete cmds;
            return -1;
        }

        if (cmds->foundOption(OPTION_BF_SIZE)) {

            bf_size = std::stoi(cmds->optionValue(OPTION_BF_SIZE));   
        }

        if (cmds->foundOption(OPTION_REQUEST_SIZE)) {

            request_size = std::stoi(cmds->optionValue(OPTION_REQUEST_SIZE));
        }

        if (cmds->foundOption(OPTION_MM_MODE)) {

            mm_mode = std::stoi(cmds->optionValue(OPTION_MM_MODE));
        }

        if (cmds->foundOption(OPTION_EH_MODE)) {

            eh_mode = std::stoi(cmds->optionValue(OPTION_EH_MODE));
        }

        if (cmds->foundOption(OPTION_RESOLV_MODE)) {

            resolv_mode = std::stoi(cmds->optionValue(OPTION_RESOLV_MODE));
        }

        if (cmds->foundOption(OPTION_ORIGIN_SERVER)) {

            strncpy(origin_server, (char *) cmds->optionValue(OPTION_ORIGIN_SERVER).c_str(), MAX_ARRAY_SIZE);
        }

        if (cmds->foundOption(OPTION_START_ROUTER)) {

            strncpy(start_router, (char *) cmds->optionValue(OPTION_START_ROUTER).c_str(), MAX_ARRAY_SIZE);
        }

        if (cmds->foundOption(OPTION_VERBOSE)) {
            verbose = true;
        }
    }

    printf("rid-analytics : request_size = %d, bf_size = %d\n", request_size, bf_size);

    if (result == ArgvParser::ParserHelpRequested) {

        delete cmds;
        return -1;
    }

    // create an rid model environment (according to the specs on nw_filename)
    RID_Analytics * rid_analytics_env = 
        new RID_Analytics(
            std::string(scn_file),          // .scn file w/ topology info
            request_size, bf_size,          // parameters for FP rate calculation
            std::string(origin_server),     // origin server location: useful for latency
            mm_mode, eh_mode, resolv_mode); // how to handle (1) multiple matches; and (2) wrong deliveries

    // ... and run the model
    rid_analytics_env->run(std::string(scn_file), std::string(start_router));

    printf("rid-analytics : output dir = %s, output label = %s\n", data_dir, output_label);

    uint8_t mode = MODE_SAVE_OUTCOMES;

    if (verbose)
        mode = (mode | MODE_VERBOSE);

    rid_analytics_env->view_results(mode, std::string(data_dir), std::string(output_label));

    // clean up after yourself...
    delete rid_analytics_env;

    return 0;
}