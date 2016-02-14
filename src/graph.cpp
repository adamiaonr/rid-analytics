#include "graph.h"

const char * GRAPHVIZ_COLORS[] = { "red", "green", "gray", "blue", "orange", "cyan", "crimson" };
const char * NODE_CHARS[] = { "MHS", "MHD", "SFP", "TPO", "DEF", "ORI", "UNK" };

Graph::Graph(const char * filename)
{
    // add .dot extension
    this->graph_filename = filename + std::string(GRAPHVIZ_DOT_FILE_EXT);

    // create the .dot output file.
    this->graph_filestream.open(this->graph_filename.c_str());

    // initial syntax to set the graph
    this->graph_filestream << "digraph {\n";
    this->graph_filestream << "\trankdir=LR\n\n";
}

Graph::~Graph() {

    this->graph_filestream.close();
}

int get_node_type_int(Node * node) {

    int node_type_int = 0;

    switch(node->get_type()) {

        case Node::MHS_NODE:
            node_type_int = MHS_NODE_INT;
            break;
        
        case Node::MHD_NODE:        
            node_type_int = MHD_NODE_INT;
            break;
        
        case Node::SFP_NODE:        
            node_type_int = SFP_NODE_INT;
            break;

        case Node::TPO_NODE:        
            node_type_int = TPO_NODE_INT;
            break;

        case Node::DEF_NODE:        
            node_type_int = DEF_NODE_INT;
            break;

        case Node::ORI_NODE:        
            node_type_int = ORI_NODE_INT;        
            break;

        default:
            node_type_int = UNK_NODE_INT;        
            break;
    }

    return node_type_int;
}

void get_node_str(char * node_str, Node * node, int depth, int breadth) {

    std::string output;

    int node_type_int = get_node_type_int(node);

    // varname
    output = NODE_CHARS[node_type_int] + std::to_string(depth) + std::to_string(3 * breadth + node_type_int);

    // FIXME: just to make it easier to add edges afterwards, use a unique 
    // id to identify each node and set it in the node directly...
    node->set_uid(output);

    // label
//    output += std::string("[label=\"") + NODE_CHARS[node_type_int] + std::to_string(depth) + std::to_string(breadth) + "\"";
//    output += std::string("[label=\"") + node->get_uid() + "\"";
    output += std::string("[label=\"") + NODE_CHARS[node_type_int] + "(" 
        + std::to_string(node->get_curr_tier()) + ":"
        + (node->get_next_tier() == END_OF_PATH ? "EOP" : std::to_string(node->get_next_tier())) 
        + ")" "\"";
    // color
    output += std::string(", color=") + GRAPHVIZ_COLORS[node_type_int];

    // wrap it up
    output += "];\n";

    snprintf(node_str, MAX_STRING_SIZE, "%s", output.c_str());
}

void get_edge_str(
    char * edge_str,
    Node * prev_node, int prev_depth, 
    Node * next_node, int next_depth,
    int breadth,
    double prob_val) {

    std::string output;

    int prev_node_type_int = get_node_type_int(prev_node);
    //int next_node_type_int = get_node_type_int(next_node);

    // initial part of edge
    // output = NODE_CHARS[get_node_type_int(prev_node)] + std::to_string(prev_depth) + std::to_string((breadth / 1)) + " -> " 
    //     + NODE_CHARS[next_node_type_int] + std::to_string(next_depth) + std::to_string((breadth * 3 + next_node_type_int)); 
    if (prev_node_type_int == ORI_NODE_INT) {

        output = NODE_CHARS[prev_node_type_int] + std::to_string(prev_depth) + std::to_string((breadth / 1)) + " -> ";

    } else {

        output = prev_node->get_uid() + " -> ";
    }
    
    output += next_node->get_uid();

    // printf("Graph::get_edge_str(): adding edge for %s(%d, %d) -> %s(%d, %d). edge_str = %s\n",
    //     NODE_CHARS[get_node_type_int(prev_node)], prev_depth, breadth,
    //     NODE_CHARS[get_node_type_int(next_node)], next_depth, (breadth * 3 + next_node_type_int),
    //     output.c_str());

    // appropriately format probability value
    char prob_val_str[MAX_STRING_SIZE];
    snprintf(prob_val_str, MAX_STRING_SIZE, "%-.3E", prob_val);
    // add probability value as 'label' and 'weight'
    output += "[label =\"" + std::string(prob_val_str) + "\", weight=\"" + std::string(prob_val_str) + "\"];\n";

    snprintf(edge_str, MAX_STRING_SIZE, "%s", output.c_str());
}

void Graph::add_node(
    Node * prev_node, int prev_depth, 
    Node * next_node, int next_depth,
    int breadth,
    double prob_val) {

    // don't add nodes which are impossible to reach...
    //if (prob_val == 0.0) return;

    char node_str[MAX_STRING_SIZE];
    char edge_str[MAX_STRING_SIZE];

    // 1) node declaration, i.e. variable name, label, color, etc.
    get_node_str(node_str, next_node, next_depth, breadth);
    this->graph_filestream << "\t" + std::string(node_str);

    // 2) edge connecting it from parent
    get_edge_str(edge_str, prev_node, prev_depth, next_node, next_depth, breadth, prob_val);
    this->graph_filestream << "\t" + std::string(edge_str) + "\n";
}

void Graph::align_nodes(
    Node ** nodes,
    int nodes_size) {

    this->graph_filestream << std::string("\t{rank = same; ");

    for (int i = 0; i < nodes_size; i++) {

        this->graph_filestream << nodes[i]->get_uid();

        if (i + 1 < nodes_size) {
            this->graph_filestream << ", ";
        }
    }    

    this->graph_filestream << " };\n"; 
}

void Graph::terminate() {

    // final bracket '{'
    this->graph_filestream << "}\n";
    this->graph_filestream.close();
}