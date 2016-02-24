#ifndef RID_ANALYTICS_HH
#define RID_ANALYTICS_HH

#define MAX_FILENAME_SIZE   64
#define MAX_REQUEST_SIZE    30
#define MAX_TIER_DEPTH      30

#define DEFAULT_ALPHA   (double) (1.0 / 10.0)
#define END_OF_PATH     (int) -1

#define AVERAGE_LATENCY     (char *) "average_latency"
#define CACHE_LATENCY       (char *) "cache_latency"
#define ORIGIN_LATENCY      (char *) "origin_latency"

#define PENALTY_TYPE_FEEDBACK   (int) 0
#define PENALTY_TYPE_FALLBACK   (int) 1

#define PENALTY_TYPE_STR_SIZE   2

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
            double * latencies;
            double * penalties;

            int * domains;          // tier breadth: "how many domains within each tier"
            int * content_sources;

            int tier_depth;
            int origin_tier;
            int cache_tier;
            int fp_resolution_tech;

            char * title;           // title for output files handled within 
                                    // (e.g. probability graphs)
        };

        RIDAnalytics() {};
        ~RIDAnalytics() {};

        static int run_model(
            double & avg_latency,
            rid_analytics_inputs input_params,
            unsigned int modes, 
            FILE ** outcomes_files,
            int outcomes_files_size,
            char * title,
            std::string data_dir);
};

#endif