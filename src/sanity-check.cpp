#include "dataparser.h"
#include "rid-analytics.h"

#define MAX_ARRAY_SIZE      128

using namespace std;

int main (int argc, char **argv) {

    double avg_latency = 0.0;
    
    // struct for input parameters
    RIDAnalytics::rid_analytics_inputs input_params;
    // these input parameters won't change during the run
    double _fp_prob[3] = {1.0E-1, 1.0E-1, 1.0E-1};
    double _alpha[3] = {0.01, 0.01, 0.01};
    double _latencies[3] = {1.0, 1.0, 1.0};
    int _domains[3] = {4, 4, 4};
    int _content_sources[3] = {0, 1, 5};

    input_params.domains = _domains;
    input_params.content_sources = _content_sources;
    input_params.fp_prob = _fp_prob;
    input_params.alpha = _alpha;
    input_params.latencies = _latencies;
    input_params.tier_depth = 3;
    input_params.fp_resolution_tech = PENALTY_TYPE_FEEDBACK;

    std::string data_dir ("test/graphs"); 

    RIDAnalytics::run_model(
        avg_latency,
        input_params,
        (MODE_VERBOSE | MODE_SAVEGRAPH),
        NULL,
        0,
        (char *) "sanity-check",
        std::string(data_dir));

    return 0;
}