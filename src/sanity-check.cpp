#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "argvparser.h"
#include "dataparser.h"
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

#define DEFAULT_SCN_FILE	(char *) "/Users/adamiaonr/workbench/rid-analytics/test/configs/sanity.scn"
#define DEFAULT_CSV_DIR     (char *) "/Users/adamiaonr/workbench/rid-analytics/test/data/sanity-check"

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
            OPTION_FWD_TABLE_SIZE,
            "size of forwarding table (in # of entries)",
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
    // keeps track of verbose mode (no verbosity by default)
    bool verbose = false;
    // evaluation parameters
    int request_size = 0;       // request size (also largest possible |F|) 
                                // in # of URL parameters
    int bf_size = 0;            // bloom filter size (in bit)
    int fwd_table_size = 0;     // forwarding table size (# of entries)


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

        if (cmds->foundOption(OPTION_FWD_TABLE_SIZE)) {

            fwd_table_size = std::stoi(cmds->optionValue(OPTION_FWD_TABLE_SIZE));
        }

        if (cmds->foundOption(OPTION_VERBOSE)) {
            verbose = true;
        }
    }

    printf("rid-analytics : request_size = %d, bf_size = %d, fwd_table_size = %d\n", 
        request_size, bf_size, fwd_table_size);

    if (result == ArgvParser::ParserHelpRequested) {

        delete cmds;
        return -1;
    }

    DataParser * scn_parser = new DataParser(scn_file);

    // if the following parameters haven't been specified via the CLI, get them 
    // from the .scn file
    if (request_size == 0)
        scn_parser->get_int_property_value(REQUEST_SIZE, request_size);

    if (bf_size == 0)
        scn_parser->get_int_property_value(BF_SIZE, bf_size);

    if (fwd_table_size == 0)
        scn_parser->get_int_property_value(FWD_TABLE_SIZE, fwd_table_size);

    // the following parameters are taken from the .scn file, only
    int access_tree_height = 0;
    scn_parser->get_int_property_value(ACCESS_TREE_HEIGHT, access_tree_height);

    int iface_num = 0;
    scn_parser->get_int_property_value(IFACE_NUM, iface_num);

    // fetch the TP list for each network router
    int tp_sizes_size = (int) pow(2, access_tree_height - 1);
    tp_sizes_size *= access_tree_height * iface_num;

    int * tp_sizes = (int *) calloc(tp_sizes_size, sizeof(int));
    int dims[3] = {access_tree_height, (int) pow(2.0, access_tree_height - 1), iface_num};

    scn_parser->get_int_property_3d_array(TP_SIZES, tp_sizes, dims);

    // fetch the iface_entry_proportion for each network router
    long double * iface_entry_proportion = 
        (long double *) calloc(
                                ((int) pow(2.0, access_tree_height - 1)) * access_tree_height * iface_num, 
                                sizeof(long double)
                            );

    scn_parser->get_double_property_3d_array(IFACE_ENTRY_PROPORTION, iface_entry_proportion, dims);

    // fetch an f_distribution and f_r_distribution set the iface info 
    // on rid_vanilla_rtr
    __float080 * f_distribution_local = (__float080 *) calloc(MAX_ARRAY_SIZE, sizeof(__float080));
    scn_parser->get_double_property_array(F_DISTRIBUTION_LOCAL, f_distribution_local);

    __float080 * f_distribution_non_local = (__float080 *) calloc(MAX_ARRAY_SIZE, sizeof(__float080));
    scn_parser->get_double_property_array(F_DISTRIBUTION_NON_LOCAL, f_distribution_non_local);

    __float080 * f_r_distribution = (__float080 *) calloc(MAX_ARRAY_SIZE, sizeof(__float080));
    scn_parser->get_double_property_array(F_R_DISTRIBUTION, f_r_distribution);    

    // FIXME: adjust f_distribution
    for (int f = 0; f < request_size; f++) {
        f_distribution_local[f] = f_distribution_local[f] / 100.0;
        f_distribution_non_local[f] = f_distribution_non_local[f] / 100.0;
        f_r_distribution[f] = f_r_distribution[f] / 100.0;    
    }

    __float080 ** f_distributions = (__float080 **) calloc(2, sizeof(__float080 *));
    f_distributions[0] = f_distribution_local;
    f_distributions[1] = f_distribution_non_local;

    RID_Analytics * rid_analytics = new RID_Analytics(
                                                1,
                                                access_tree_height,
                                                iface_num,
                                                request_size,
                                                bf_size,
                                                fwd_table_size,
                                                iface_entry_proportion,
                                                f_distributions);

    rid_analytics->run(request_size, tp_sizes, f_r_distribution);

    char output_file_path[MAX_ARRAY_SIZE] = {0};
    snprintf(output_file_path, MAX_ARRAY_SIZE, "%s/%s.csv", data_dir, output_file);

    printf("rid-analytics : output_file_path = %s\n", output_file_path);

    uint8_t mode = MODE_SAVE_OUTCOMES;

    if (verbose)
        mode = (mode | MODE_VERBOSE);

    rid_analytics->view_results(mode, output_file_path);

    // clean up after yourself...
    delete scn_parser;
    delete rid_analytics;

    free(tp_sizes);
    free(iface_entry_proportion);
    free(f_distribution_local);
    free(f_distribution_non_local);
    free(f_distributions);
    free(f_r_distribution);

    return 0;
}