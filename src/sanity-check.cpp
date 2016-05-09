#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "dataparser.h"
#include "rid-analytics.h"

#define MAX_ARRAY_SIZE      128

#define REQUEST_SIZE                (char *) "request_size"
#define IFACE_ENTRY_PROPORTION      (char *) "iface_entry_proportion_"
#define F_DISTRIBUTION              (char *) "f_distribution"
#define F_R_DISTRIBUTION            (char *) "f_r_distribution"
#define ACCESS_TREE_HEIGHT          (char *) "access_tree_height"
#define TP_SIZES                    (char *) "tp_size_"
#define IFACE_NUM                   (char *) "iface_num"
#define FWD_TABLE_SIZE              (char *) "fwd_table_size"

#define DEFAULT_SCN_FILE	(char *) "/Users/adamiaonr/workbench/rid-analytics/test/configs/sanity.scn"
#define DEFAULT_CSV_DIR     (char *) "/Users/adamiaonr/workbench/rid-analytics/test/data/sanity-check"

using namespace std;

int main (int argc, char **argv) {

    DataParser * scn_parser = new DataParser(DEFAULT_SCN_FILE);

    // fetch access network parameters
    int request_size = 0;
    scn_parser->get_int_property_value(REQUEST_SIZE, request_size);

    int access_tree_height = 0;
    scn_parser->get_int_property_value(ACCESS_TREE_HEIGHT, access_tree_height);

    int iface_num = 0;
    scn_parser->get_int_property_value(IFACE_NUM, iface_num);

    int fwd_table_size = 0;
    scn_parser->get_int_property_value(FWD_TABLE_SIZE, fwd_table_size);

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
    __float080 * f_distribution = (__float080 *) calloc(MAX_ARRAY_SIZE, sizeof(__float080));
    scn_parser->get_double_property_array(F_DISTRIBUTION, f_distribution);

    __float080 * f_r_distribution = (__float080 *) calloc(MAX_ARRAY_SIZE, sizeof(__float080));
    scn_parser->get_double_property_array(F_R_DISTRIBUTION, f_r_distribution);    

    // FIXME: adjust f_distribution
    for (int f = 0; f < request_size; f++) {
        f_distribution[f] = f_distribution[f] / 100.0;
        f_r_distribution[f] = f_r_distribution[f] / 100.0;    
    }

    RID_Analytics * rid_analytics = new RID_Analytics(
                                                1,
                                                access_tree_height,
                                                iface_num,
                                                request_size,
                                                fwd_table_size,
                                                iface_entry_proportion,
                                                f_distribution);

    rid_analytics->run(request_size, tp_sizes, f_r_distribution);

    char output_file_path[MAX_ARRAY_SIZE] = {0};
    snprintf(output_file_path, MAX_ARRAY_SIZE, "%s/no0-4h.csv", DEFAULT_CSV_DIR);

    rid_analytics->view_results(MODE_VERBOSE | MODE_SAVE_OUTCOMES, output_file_path);

    // clean up after yourself...
    delete scn_parser;
    delete rid_analytics;

    free(tp_sizes);
    free(iface_entry_proportion);
    free(f_distribution);
    free(f_r_distribution);

    return 0;
}