#ifndef NODE_HH
#define NODE_HH

#define MAX_NODE_STRING_SIZE 128

#define OUTCOME_CCACHE              (char *) "correct cache"
#define OUTCOME_IDEST_CSERVER       (char *) "wrong dest > server"
#define OUTCOME_DROPPED             (char *) "dropped > server"

class Node {

    public:

        enum Type {
            C_NODE = 0x00, 
            I_NODE = 0x01,
            N_NODE = 0x02,
            O_NODE = 0x04,
            UNKNOWN = 0x08
        };

        Node() {

            this->curr_level = 0;
            this->next_level = 0;
            this->max_level = 0;
            this->prob_val = 0.0;
            this->latency_val = 0.0;
            this->type = Node::UNKNOWN;
        }

        Node(Node::Type type) {

            this->curr_level = 0;
            this->next_level = 0;
            this->max_level = 0;
            this->prob_val = 0.0;
            this->latency_val = 0.0;
            this->type = type;
        }

        Node(int max_level, int curr_level, Node::Type type) {

            this->curr_level = curr_level;
            this->next_level = 0;
            this->max_level = max_level;
            this->prob_val = 0.0;
            this->latency_val = 0.0;
            this->type = type;
        }

        Node(int curr_level, int next_level, int max_level, double prob_val, double latency_val) {

            this->curr_level = curr_level;
            this->next_level = next_level;
            this->max_level = max_level;
            this->prob_val = prob_val;
            this->latency_val = latency_val;
            this->type = Node::UNKNOWN;
        }

        Node(int curr_level, int next_level, int max_level, double prob_val, double latency_val, Node::Type type) {

            this->curr_level = curr_level;
            this->next_level = next_level;
            this->max_level = max_level;
            this->prob_val = prob_val;
            this->latency_val = latency_val;
            this->type = type;
        }

        ~Node() {};

        void set_curr_level(int curr_level) { this->curr_level = curr_level; }
        void set_next_level(int next_level) { this->next_level = next_level; }
        void set_max_level(int max_level) { this->max_level = max_level; }
        void set_prob_val(double prob_val) { this->prob_val = prob_val; }
        void set_latency_val(double latency_val) { this->latency_val = latency_val; }
        void set_type(Node::Type type) { this->type = type; }
        void set_uid(std::string uid) { this->uid = uid; }
        void set_outcome(std::string outcome) { this->outcome = outcome; }

        int get_curr_level() const { return this->curr_level; }
        int get_next_level() const { return this->next_level; }
        int get_max_level() const { return this->max_level; }

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

                case Node::I_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "I");
                    break;
                case Node::C_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "C");
                    break;
                case Node::N_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "N");
                    break;
                case Node::O_NODE:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "O");
                    break;
                default:
                    snprintf(type_str, MAX_NODE_STRING_SIZE, "U");
                    break;
            }

            return type_str;
        }

        char * to_string() {

            // using calloc() since node_str will be returned
            char * node_str = (char *) calloc(MAX_NODE_STRING_SIZE, sizeof(char));
            char node_type_str[MAX_NODE_STRING_SIZE];

            switch(this->type) {

                case Node::I_NODE:
                    snprintf(node_type_str, MAX_NODE_STRING_SIZE, "I_NODE");
                    break;
                case Node::C_NODE:
                    snprintf(node_type_str, MAX_NODE_STRING_SIZE, "C_NODE");
                    break;
                case Node::N_NODE:
                    snprintf(node_type_str, MAX_NODE_STRING_SIZE, "N_NODE");
                    break;
                case Node::O_NODE:
                    snprintf(node_type_str, MAX_NODE_STRING_SIZE, "O_NODE");
                    break;
                default:
                    snprintf(node_type_str, MAX_NODE_STRING_SIZE, "O_UNKNOWN");
                    break;
            }

            snprintf(
                node_str, MAX_NODE_STRING_SIZE, 
                "[TYPE : %s ][PROB : %02.6E][LATENCY : %08.6f][MAX_LEVEL : %d][OUTCOME: %s]",
                node_type_str, this->prob_val, this->latency_val, this->max_level, this->outcome.c_str());

            return node_str;
        }

    private:

        int curr_level;
        int next_level;
        int max_level;
        double prob_val;
        double latency_val;
        Node::Type type;
        std::string uid;
        std::string outcome;
};

#endif