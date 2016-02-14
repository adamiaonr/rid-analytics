#include "dataparser.h"
#include "rid-analytics.h"

#define MAX_ARRAY_SIZE      128

using namespace std;

int main (int argc, char **argv) {

    double avg_latency = 0.0;
    
    // struct for input parameters
    RIDAnalytics::rid_analytics_inputs input_params;
    // these input parameters won't change during the run
    double _fp_prob[4] = {0.001, 0.001, 1.0, 0.001};
    double _alpha[4] = {0.5, 0.5, 0.5, 0.5};
    int _domains[4] = {4, 4, 4, 4};
    int _content_sources[4] = {0, 1, 1, 1};

    input_params.domains = _domains;
    input_params.content_sources = _content_sources;
    input_params.fp_prob = _fp_prob;
    input_params.alpha = _alpha;
    input_params.tier_depth = 4;
    input_params.cache_tier = 1;
    input_params.origin_tier = 3;
    input_params.fp_resolution_tech = PENALTY_TYPE_FEEDBACK;

    std::string data_dir ("test/graphs"); 

    RIDAnalytics::run_model(
        avg_latency,
        input_params,
        (MODE_VERBOSE | MODE_SAVEGRAPH),
        NULL,
        std::string(data_dir));

    return 0;
}