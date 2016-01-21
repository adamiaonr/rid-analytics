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

#define I_NODE_INT 0x00
#define C_NODE_INT 0x01
#define N_NODE_INT 0x02
#define O_NODE_INT 0x03
#define U_NODE_INT 0x04

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
            Node * i_node,
            Node * c_node,
            Node * n_node);

    private: 

        std::string graph_filename;
        ofstream graph_filestream;
};

#endif