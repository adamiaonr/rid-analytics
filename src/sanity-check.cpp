#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "argvparser.h"
#include "rid-analytics.h"

#define MAX_ARRAY_SIZE      256

#define OPTION_SCN_FILE             (char *) "scn-file"
#define OPTION_DATA_DIR             (char *) "data-dir"
#define OPTION_OUTPUT_FILE          (char *) "output-file"
#define OPTION_VERBOSE              (char *) "verbose"
// eval parameters can be specified directly through the CLI and overwrite any 
// parameters set in .scn files
#define OPTION_REQUEST_SIZE         (char *) "request-size"
#define OPTION_BF_SIZE              (char *) "bf-size"
#define OPTION_FWD_TABLE_SIZE       (char *) "fwd-table-size"
#define OPTION_F_MIN_ANNC           (char *) "f-min-annc"
// special eval options : 
//  * --f-min : the min. length of forwarding entries allowed in the forwarding 
//              table
//  * --expand-factor : expanding entries of length f to (f + a) creates 
//                      f * (expand-factor)^(a) new entries
#define OPTION_F_MIN                (char *) "f-min"
#define OPTION_EXPAND_FACTOR        (char *) "expand-factor"

#define REQUEST_SIZE                (char *) "request_size"
#define BF_SIZE                     (char *) "bf_size"
#define IFACE_ENTRY_PROPORTION      (char *) "iface_entry_proportion_"
#define F_DISTRIBUTION_LOCAL        (char *) "f_distribution_local"
#define F_DISTRIBUTION_NON_LOCAL    (char *) "f_distribution_non_local"
#define F_R_DISTRIBUTION            (char *) "f_r_distribution"
#define ACCESS_TREE_HEIGHT          (char *) "access_tree_height"
#define TP_SIZES                    (char *) "tp_size_"
#define IFACE_NUM                   (char *) "iface_num"
#define FWD_TABLE_SIZE              (char *) "fwd_table_size"

#define NON_LOCAL_FRACTION          (__float080) 0.75
#define LOCAL_FRACTION              (__float080) 0.25

using namespace std;
using namespace CommandLineProcessing;

ArgvParser * create_argv_parser() {

    ArgvParser * cmds = new ArgvParser();

    cmds->setIntroductoryDescription("\n\nrid-analytics v1.0\n\na tool to simulate the behavior "\
        "of RID networks. it computes the probability of different delivery states for "\
        "all possible request paths in a network topology passed as input.\nby adamiaonr@cmu.edu");

    cmds->setHelpOption("h", "help", "help page.");

    cmds->defineOption(
            OPTION_SCN_FILE,
            "path to .scn file which contains the input info for the analysis",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_DATA_DIR,
            "path to a directory to output .csv files w/ latency data",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_OUTPUT_FILE,
            "name of .csv file to gather analysis data (will be created in <data-dir>)",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_BF_SIZE,
            "size of RID or Bloom filter (in bits)",
            ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_REQUEST_SIZE,
            "size of request, in # of URL components",
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
    char output_file[MAX_ARRAY_SIZE];
    // verbose mode (no verbosity by default)
    bool verbose = false;
    // evaluation parameters
    // request size (also largest possible |F|) in # of URL parameters
    int request_size = 0;       
    // bloom filter size (in bit)
    int bf_size = 0;

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

        if (cmds->foundOption(OPTION_OUTPUT_FILE)) {

            strncpy(output_file, (char *) cmds->optionValue(OPTION_OUTPUT_FILE).c_str(), MAX_ARRAY_SIZE);

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
    RID_Analytics * rid_analytics_env = new RID_Analytics(std::string(scn_file), request_size, bf_size);
    // ... and run the model
    rid_analytics_env->run(std::string(scn_file));

    printf("rid-analytics : output dir = %s\n", data_dir);

    uint8_t mode = 0x00;
    if (verbose)
        mode = (mode | MODE_VERBOSE);

    rid_analytics_env->view_results(mode, std::string(data_dir));

    // clean up after yourself...
    delete rid_analytics_env;

    return 0;
}