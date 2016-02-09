#ifndef RID_ANALYTICS_HH
#define RID_ANALYTICS_HH

#define MAX_FILENAME_SIZE   64
#define MAX_REQUEST_SIZE    30
#define MAX_TIER_DEPTH     30

#define DEFAULT_O_OPTIMISTIC (double) (1.0 / 10.0)
#define END_OF_PATH (int) -1

#define AVERAGE_LATENCY     (char *) "average_latency"
#define CACHE_LATENCY       (char *) "cache_latency"
#define ORIGIN_LATENCY      (char *) "origin_latency"

#define PENALTY_TYPE_FEEDBACK   0x00
#define PENALTY_TYPE_FALLBACK   0x01

#define MODE_VERBOSE        0x01
#define MODE_SAVECDF        0x02
#define MODE_SAVEGRAPH      0x04
#define MODE_SAVEOUTCOMES   0x08

#define BF_SIZE 160

class RIDAnalytics {

    public:

        enum DecisionMode {
            DEFAULT = 0x00
        };

        struct rid_analytics_inputs {

            double * fp_prob;
            double * alpha;

            int tier_depth;
            int cache_tier;
            int origin_tier;
            int fp_resolution_tech;

        };

        RIDAnalytics() {};
        ~RIDAnalytics() {};

        static int run_model(
            double & avg_latency,
            rid_analytics_inputs input_params,
            unsigned int modes, 
            FILE ** outcomes_file,
            std::string data_dir);
};

#endif