#ifndef RID_ANALYTICS_HH
#define RID_ANALYTICS_HH

#define MAX_FILENAME_SIZE   64
#define MAX_REQUEST_SIZE    30
#define MAX_LEVEL_DEPTH     30

#define DEFAULT_O_OPTIMISTIC (double) (1.0 / 10.0)
#define END_OF_PATH (int) -1

#define AVERAGE_LATENCY     (char *) "average_latency"
#define CACHE_LATENCY       (char *) "cache_latency"
#define ORIGIN_LATENCY      (char *) "origin_latency"

#define PENALTY_TYPE_FEEDBACK   0x00
#define PENALTY_TYPE_FALLBACK   0x01

#define FP_RATE_TYPE_PRMT           0x00
#define FP_RATE_TYPE_AUTO_V         0x01
#define FP_RATE_TYPE_AUTO_INV_V     0x02
#define FP_RATE_TYPE_AUTO_RISE      0x03
#define FP_RATE_TYPE_AUTO_DECREASE  0x04
#define FP_RATE_TYPE_AUTO_HIGH      0x05
#define FP_RATE_TYPE_AUTO_LOW       0x06

#define BF_SIZE 160

class RIDAnalytics {

    public:

        enum DecisionMode {
            DEFAULT = 0x00
        };

        RIDAnalytics() {};
        ~RIDAnalytics() {};

        static int run_model(
            double & avg_latency,
            double * fp_prob, 
            double * o_optimistic, 
            int level_depth,
            int cache_level, 
            int origin_level,
            int fp_resolution_tech,
            bool save_cdf,
            bool save_graph,
            std::string data_dir);
};

#endif