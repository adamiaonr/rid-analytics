#include "dataparser.h"
#include "rid-analytics.h"

#define MAX_ARRAY_SIZE      128

using namespace std;

int main (int argc, char **argv) {

    double avg_latency = 0.0;
    double _fp_prob[4] = {0.001, 0.001, 1.0, 0.001};
    double _o_optimistic[4] = {1.0, 1.0, 1.0, 1.0};
    
    int levels = 4;
    int c_level = 1;
    int o_level = 3;

    std::string data_dir ("test/graphs"); 

    RIDAnalytics::run_model(
        avg_latency,
        _fp_prob, 
        _o_optimistic, 
        levels,
        c_level, 
        o_level,
        PENALTY_TYPE_FEEDBACK,
        true,
        false,
        true,
        false,
        NULL,
        NULL,
        data_dir);

    return 0;
}