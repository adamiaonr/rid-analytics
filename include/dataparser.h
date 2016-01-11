#ifndef DATA_PARSER_HH
#define DATA_PARSER_HH

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

#define MAX_FILENAME_SIZE 64
#define MAX_CHARS_PER_LINE 512
#define MAX_TOKENS_PER_LINE 20

#define PROPERTY_DELIMITER (const char *) "="
#define VALUE_DELIMITER (const char *) ","

using namespace std;

class DataParser {

    public:

        DataParser() {}
        DataParser(const char * dat_filename);
        ~DataParser() {}

        int get_int_property_value(const char * property_name, int &int_value);
        int get_int_property_array(const char * property_name, int * int_array);

        int get_double_property_value(const char * property_name, double &double_value);
        int get_double_property_array(const char * property_name, double * double_array);

        int get_string_property_value(const char * property_name, std::string &string_value);
        int get_string_property_array(const char * property_name, std::string * string_array);

    private: 

        std::string dat_filename;
};

#endif