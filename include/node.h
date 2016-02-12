#ifndef NODE_HH
#define NODE_HH

#define MAX_NODE_STRING_SIZE 128

#define OUTCOME_CCACHE              (char *) "correct cache"
#define OUTCOME_IDEST_CSERVER       (char *) "wrong dest > server"
#define OUTCOME_DROPPED             (char *) "dropped > relay"

class Node {

    public:

        enum Type {
            MHS_NODE = 0x00, 
            MHD_NODE = 0x01,
            SFP_NODE = 0x02,
            TPO_NODE = 0x04,
            DEF_NODE = 0x08,
            ORI_NODE = 0x10,
            UNKNOWN = 0x20
        };

        Node() {

            this->curr_tier = 0;
            this->next_tier = 0;
            this->max_tier = 0;
            this->prob_val = 0.0;
            this->latency_val = 0.0;
            this->type = Node::UNKNOWN;
        }

        Node(Node::Type type) {

            this->curr_tier = 0;
            this->next_tier = 0;
            this->max_tier = 0;
            this->prob_val = 0.0;
            this->latency_val = 0.0;
            this->type = type;
        }

        Node(int max_tier, int curr_tier, Node::Type type) {

            this->curr_tier = curr_tier;
            this->next_tier = 0;
            this->max_tier = max_tier;
            this->prob_val = 0.0;
            this->latency_val = 0.0;
            this->type = type;
        }

        Node(int curr_tier, int next_tier, int max_tier, double prob_val, double latency_val) {

            this->curr_tier = curr_tier;
            this->next_tier = next_tier;
            this->max_tier = max_tier;
            this->prob_val = prob_val;
            this->latency_val = latency_val;
            this->type = Node::UNKNOWN;
        }

        Node(int curr_tier, int next_tier, int max_tier, double prob_val, double latency_val, Node::Type type) {

            this->curr_tier = curr_tier;
            this->next_tier = next_tier;
            this->max_tier = max_tier;
            this->prob_val = prob_val;
            this->latency_val = latency_val;
            this->type = type;
        }

        ~Node() {};

        void set_curr_tier(int curr_tier) { this->curr_tier = curr_tier; }
        void set_next_tier(int next_tier) { this->next_tier = next_tier; }
        void set_max_tier(int max_tier) { this->max_tier = max_tier; }
        void set_prob_val(double prob_val) { this->prob_val = prob_val; }
        void set_latency_val(double latency_val) { this->latency_val = latency_val; }
        void set_type(Node::Type type) { this->type = type; }
        void set_uid(std::string uid) { this->uid = uid; }
        void set_outcome(std::string outcome) { this->outcome = outcome; }

        int get_curr_tier() const { return this->curr_tier; }
        int get_next_tier() const { return this->next_tier; }
        int get_max_tier() const { return this->max_tier; }

        double get_prob_val() const { return this->prob_val; }
        double get_latency_val() const { return this->latency_val; }

        // check this stackoverflow post for a 'const' explanation: 
        // http://goo.gl/0u4Alx
        // TL;DR: gcc assumes that member functions *NOT* declared as 
        // const *WILL* modify const objects, and so terminates compilation with 
        // an error. solution: make the functions const.
        std::string get_uid() const { return this->uid; }
        std::string get_outcome() const { return this->outcome; }

        Node::Type get_type() const { return this->type; }

        char * get_type_str() {
        
            char * type_str = (char *) calloc(MAX_NODE_STRING_SIZE, sizeof(char));

            switch(this->type) {

                case Node::MHS_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "MHS");
                    break;
                case Node::MHD_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "MHD");
                    break;
                case Node::SFP_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "SFP");
                    break;
                case Node::TPO_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "TPO");
                    break;
                case Node::DEF_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "DEF");
                    break;
                case Node::ORI_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "ORI");
                    break;
                default:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "UNK");
                    break;
            }

            return type_str;
        }

        char * to_string() {

            // using calloc() since node_str will be returned
            char * node_str = (char *) calloc(MAX_NODE_STRING_SIZE, sizeof(char));
            char type_str[MAX_NODE_STRING_SIZE];

            switch(this->type) {

                case Node::MHS_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "MHS");
                    break;
                case Node::MHD_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "MHD");
                    break;
                case Node::SFP_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "SFP");
                    break;
                case Node::TPO_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "TPO");
                    break;
                case Node::DEF_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "DEF");
                    break;
                case Node::ORI_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "ORI");
                    break;
                default:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "UNK");
                    break;
            }

            snprintf(
                node_str, MAX_NODE_STRING_SIZE, 
                "[TYPE : %s ][PROB : %02.6E][LATENCY : %08.6f][MAX_TIER : %d][OUTCOME: %s]",
                type_str, this->prob_val, this->latency_val, this->max_tier, this->outcome.c_str());

            return node_str;
        }

    private:

        int curr_tier;
        int next_tier;
        int max_tier;
        double prob_val;
        double latency_val;
        Node::Type type;
        std::string uid;
        std::string outcome;
};

#endif