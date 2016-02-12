#ifndef GRAPH_HH
#define GRAPH_HH

#include <dirent.h>
#include <errno.h>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <vector>
#include <math.h>
#include <cmath>
#include <cstring>
#include <sys/time.h>

#include "node.h"

#define MAX_STRING_SIZE 128

#define GRAPHVIZ_DOT_FILE_EXT (const char *) ".dot"

#define MHS_NODE_INT 0x00
#define MHD_NODE_INT 0x01
#define SFP_NODE_INT 0x02
#define TPO_NODE_INT 0x04
#define DEF_NODE_INT 0x08
#define ORI_NODE_INT 0x08
#define UNKNOWN_INT 0x08

#define END_OF_PATH (int) -1

using namespace std;

class Graph {

    public:

        Graph() {}
        Graph(const char * dat_filename);
        ~Graph();

        void terminate();
        
        void add_node(
            Node * prev_node, int prev_level, 
            Node * next_node, int next_level,
            int breadth,
            double prob_val);
        
        void align_nodes(
            Node ** nodes,
            int nodes_size);

    private: 

        std::string graph_filename;
        ofstream graph_filestream;
};

#endif