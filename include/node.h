#ifndef NODE_HH
#define NODE_HH

#define MAX_NODE_STRING_SIZE 128

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

            this->level = 0;
            this->max_level = 0;
            this->prob_val = 0.0;
            this->latency_val = 0.0;
            this->type = Node::UNKNOWN;
        }

        Node(Node::Type type) {

            this->level = 0;
            this->max_level = 0;
            this->prob_val = 0.0;
            this->latency_val = 0.0;
            this->type = type;
        }

        Node(int max_level, Node::Type type) {

            this->level = 0;
            this->max_level = max_level;
            this->prob_val = 0.0;
            this->latency_val = 0.0;
            this->type = type;
        }

        Node(int level, int max_level, double prob_val, double latency_val) {

            this->level = level;
            this->max_level = max_level;
            this->prob_val = prob_val;
            this->latency_val = latency_val;
            this->type = Node::UNKNOWN;
        }

        Node(int level, int max_level, double prob_val, double latency_val, Node::Type type) {

            this->level = level;
            this->max_level = max_level;
            this->prob_val = prob_val;
            this->latency_val = latency_val;
            this->type = type;
        }

        ~Node() {};

        void set_level(int level) { this->level = level; }
        void set_max_level(int max_level) { this->max_level = max_level; }
        void set_type(Node::Type type) { this->type = type; }
        void set_prob_val(double prob_val) { this->prob_val = prob_val; }
        void set_latency_val(double latency_val) { this->latency_val = latency_val; }
        void set_uid(std::string uid) { this->uid = uid; }

        int get_level() { return this->level; }
        int get_max_level() { return this->max_level; }
        double get_prob_val() { return this->prob_val; }
        double get_latency_val() { return this->latency_val; }

        std::string get_uid() { return this->uid; }

        Node::Type get_type() { return this->type; }

        char * to_string() {

            // using calloc() since node_str will be returned
            char * node_str = (char *) (char *) calloc(MAX_NODE_STRING_SIZE, sizeof(char));
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
                "[TYPE : %s ][PROB : %02.6E][LATENCY : %08.6f][MAX_LEVEL : %d]",
                node_type_str, this->prob_val, this->latency_val, this->max_level);

            return node_str;
        }

    private:

        int level;
        int max_level;
        double prob_val;
        double latency_val;
        Node::Type type;

        std::string uid;
};

#endif