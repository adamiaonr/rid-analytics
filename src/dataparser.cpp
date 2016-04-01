#include "dataparser.h"

#include <iostream>

DataParser::DataParser(const char * dat_filename)
{
    this->dat_filename = std::string(dat_filename);
}

int DataParser::get_int_property_value(const char * property_name, int &int_value) 
{
    std::ifstream ifs;
    ifs.open(this->dat_filename.c_str());

    if (!ifs.good()) {

        fprintf(
            stderr, 
            "DataParser::get_int_property_value() : .scn file (%s) not found.\n", 
            this->dat_filename.c_str());

        ifs.close();

        return -1;
    }

    while (!ifs.eof()) {

        // read line of the .scn file
        char buf[MAX_CHARS_PER_LINE];
        ifs.getline(buf, MAX_CHARS_PER_LINE);

        // array to store memory addresses of the tokens in buf
        const char *token[MAX_TOKENS_PER_LINE] = {}; 

        // extract the first token, the property name
        token[0] = strtok(buf, PROPERTY_DELIMITER);

        if (strcmp(token[0], property_name) != 0) {
            // ...this isn't the property we're looking for. move along!
            continue;
        }

        // since this is a single int token, simply look for the next token and 
        // 'atoi()' it
        if (token[0]) {

            token[1] = strtok(0, VALUE_DELIMITER);

            if (token[1]) {

                int_value = atoi(token[1]);
                ifs.close();

                return 0;

            } else {

                fprintf(
                    stderr, 
                    "DataParser::get_int_property_value() : no value for property %s\n", 
                    token[0]);

                ifs.close();

                return -1;
            }
        }
    }

    fprintf(
        stderr, 
        "DataParser::get_int_property_value() : no property w/ name %s\n", 
        property_name);

    ifs.close();

    return -1;
}

int DataParser::get_int_property_array(const char * property_name, int * int_array) 
{
    std::ifstream ifs;
    ifs.open(this->dat_filename.c_str());

    if (!ifs.good()) {

        fprintf(
            stderr, 
            "DataParser::get_int_property_array() : .scn file (%s) not found.\n", 
            this->dat_filename.c_str());

        ifs.close();

        return -1;
    }

    int n = 0;

    while (!ifs.eof()) {

        // read line of the .scn file
        char buf[MAX_CHARS_PER_LINE];
        ifs.getline(buf, MAX_CHARS_PER_LINE);

        // array to store memory addresses of the tokens in buf
        const char *token[MAX_TOKENS_PER_LINE] = {}; 

        // extract the first token, the property name
        token[0] = strtok(buf, PROPERTY_DELIMITER);

        if (strcmp(token[0], property_name) != 0) {
            // ...this isn't the property we're looking for. move along!
            continue;
        }

        //fprintf(stderr, "DataParser::get_int_property_array() : got property = %s\n", token[0]);

        if (token[0]) {

            for (n = 0; n < MAX_TOKENS_PER_LINE; ++n) {

                token[n] = strtok(0, VALUE_DELIMITER);

                if (!token[n]) {

                    break;

                } else {

                    //fprintf(stderr, "DataParser::get_int_property_array() : got value = %s\n", token[n]);

                    int_array[n] = atoi(token[n]);
                }
            }

            ifs.close();

            if (n > 0) {

                return 0;

            } else {

                fprintf(
                    stderr, 
                    "DataParser::get_int_property_array() : no tokens gathered\n");

                return -1;
            }
        }
    }

    ifs.close();

    fprintf(
        stderr, 
        "DataParser::get_int_property_value() : no property w/ name %s\n", property_name);

    return -1;
}

int DataParser::get_int_property_3d_array(
    const char * property_prefix, 
    int * int_array, 
    int * dim) 
{
    // temporarily fetch the whole row of index x of the 3D array
    if (int_array == NULL) {

        fprintf(
            stderr, 
            "DataParser::get_int_property_3d_array() : int_array is NULL.\n");

        return -1;
    }

    int _int_array[1024] = {0};

    char property_name[128];
    char x_str[8];

    uint32_t _offset_x = 0;
    uint32_t _offset_y = 0;

    for (uint8_t x = 0; x < dim[0]; x++) {

        // clean the row array
        memset(_int_array, 0, 1024);

        // clean up property_name and x_str
        memset(property_name, 0, 128);
        memset(x_str, 0, 8);

        // re-construct property_name : start with property_prefix
        strncpy(property_name, property_prefix, 128);

        // concat the property suffix
        snprintf(x_str, 8, "%d", x);
        strncat(property_name, x_str, 128);

        // temporarily fetch the whole row of index x of the 3D array
        if (get_int_property_array(property_name, _int_array) < 0) {

            fprintf(
                stderr, 
                "DataParser::get_int_property_3d_array() : error while reading row %s.\n", 
                property_name);

            return -1;
        }

        _offset_x = x * (dim[1] * dim[2]);

        for (uint8_t y = 0; y < dim[1]; y++) {

            _offset_y = y * dim[2];

            for (uint8_t z = 0; z < dim[2]; z++)
                int_array[_offset_x + _offset_y + z] = _int_array[_offset_y + z];
        }
    }

    return 0;
}

int DataParser::get_double_property_value(const char * property_name, long double &double_value) 
{
    std::ifstream ifs;
    ifs.open(this->dat_filename.c_str());

    if (!ifs.good()) {

        fprintf(
            stderr, 
            "DataParser::get_double_property_value() : .scn file (%s) not found.\n", 
            this->dat_filename.c_str());

        ifs.close();

        return -1;
    }

    while (!ifs.eof()) {

        // read line of the .scn file
        char buf[MAX_CHARS_PER_LINE];
        ifs.getline(buf, MAX_CHARS_PER_LINE);

        // array to store memory addresses of the tokens in buf
        const char *token[MAX_TOKENS_PER_LINE] = {}; 

        // extract the first token, the property name
        token[0] = strtok(buf, PROPERTY_DELIMITER);

        if (strcmp(token[0], property_name) != 0) {
            // ...this isn't the property we're looking for. move along!
            continue;
        }

        // since this is a single int token, simply look for the next token and 
        // 'atoi()' it
        if (token[0]) {

            token[1] = strtok(0, VALUE_DELIMITER);

            if (token[1]) {

                double_value = atof(token[1]);
                ifs.close();

                return 0;

            } else {

                fprintf(
                    stderr, 
                    "DataParser::get_double_property_value() : no value for property %s\n", 
                    token[0]);

                ifs.close();

                return -1;
            }
        }
    }

    fprintf(
        stderr, 
        "DataParser::get_double_property_value() : no property w/ name %s\n", 
        property_name);

    ifs.close();

    return -1;
}

int DataParser::get_double_property_array(
    const char * property_name, 
    long double * double_array) 
{
    std::ifstream ifs;
    ifs.open(this->dat_filename.c_str());

    if (!ifs.good()) {

        fprintf(
            stderr, 
            "DataParser::get_double_property_array() : .scn file (%s) not found.\n", 
            this->dat_filename.c_str());

        ifs.close();

        return -1;
    }

    int n = 0;

    while (!ifs.eof()) {

        // read line of the .scn file
        char buf[MAX_CHARS_PER_LINE];
        ifs.getline(buf, MAX_CHARS_PER_LINE);

        // array to store memory addresses of the tokens in buf
        const char *token[MAX_TOKENS_PER_LINE] = {}; 

        // extract the first token, the property name
        token[0] = strtok(buf, PROPERTY_DELIMITER);

        if (token[0]) {

            if (strcmp(token[0], property_name) != 0) {
                // ...this isn't the property we're looking for. move along!
                continue;
            }

            for (n = 0; n < MAX_TOKENS_PER_LINE; ++n) {

                token[n] = strtok(0, VALUE_DELIMITER);

                if (!token[n]) {

                    break;

                } else {

                    double_array[n] = (long double) atof(token[n]);
                }
            }

            ifs.close();

            if (n > 0) {

                return 0;

            } else {

                fprintf(
                    stderr, 
                    "DataParser::get_double_property_array() : no tokens gathered\n");

                return -1;
            }
        }
    }

    ifs.close();

    fprintf(
        stderr, 
        "DataParser::get_double_property_array() : no property w/ name %s\n", 
        property_name);

    return -1;
}

int DataParser::get_double_property_3d_array(
    const char * property_prefix, 
    long double * double_array, 
    int * dim) 
{

    if (double_array == NULL) {

        fprintf(
            stderr, 
            "DataParser::get_double_property_3d_array() : double_array is NULL.\n");

        return -1;
    }

    long double _double_array[2048] = {0};

    char property_name[128];
    char x_str[8];

    uint32_t _offset_x = 0;
    uint32_t _offset_y = 0;

    for (uint8_t x = 0; x < dim[0]; x++) {

        // clean the row array
        memset(_double_array, 0, 1024);

        // clean up property_name and x_str
        memset(property_name, 0, 128);
        memset(x_str, 0, 8);

        // re-construct property_name : start with property_prefix
        strncpy(property_name, property_prefix, 128);

        // concat the property suffix
        snprintf(x_str, 8, "%d", x);
        strncat(property_name, x_str, 128);

        // temporarily fetch the whole row of index x of the 3D array
        if (get_double_property_array(property_name, _double_array) < 0) {

            fprintf(
                stderr, 
                "DataParser::get_double_property_3d_array() : error while reading row %s.\n", 
                property_name);

            return -1;
        }

        _offset_x = x * (dim[1] * dim[2]);

        for (uint8_t y = 0; y < dim[1]; y++) {

            _offset_y = y * dim[2];

            for (uint8_t z = 0; z < dim[2]; z++) {
                double_array[_offset_x + _offset_y + z] = _double_array[_offset_y + z];
            }
        }
    }

    return 0;
}

int DataParser::get_string_property_value(
    const char * property_name, 
    std::string &string_value) 
{
    std::ifstream ifs;
    ifs.open(this->dat_filename.c_str());

    if (!ifs.good()) {

        fprintf(
            stderr, 
            "DataParser::get_string_property_value() : .scn file (%s) not found.\n", 
            this->dat_filename.c_str());

        ifs.close();

        return -1;
    }

    while (!ifs.eof()) {

        // read line of the .scn file
        char buf[MAX_CHARS_PER_LINE];
        ifs.getline(buf, MAX_CHARS_PER_LINE);

        // array to store memory addresses of the tokens in buf
        const char *token[MAX_TOKENS_PER_LINE] = {}; 

        // extract the first token, the property name
        token[0] = strtok(buf, PROPERTY_DELIMITER);

        if (strcmp(token[0], property_name) != 0) {
            // ...this isn't the property we're looking for. move along!
            continue;
        }

        // since this is a single int token, simply look for the next token and 
        // 'atoi()' it
        if (token[0]) {

            token[1] = strtok(0, VALUE_DELIMITER);

            if (token[1]) {

                string_value = std::string(token[1]);
                ifs.close();

                return 0;

            } else {

                fprintf(
                    stderr, 
                    "DataParser::get_string_property_value() : no value for property %s\n", 
                    token[0]);

                ifs.close();

                return -1;
            }
        }
    }

    fprintf(
        stderr, 
        "DataParser::get_string_property_value() : no property w/ name %s\n", 
        property_name);

    ifs.close();

    return -1;
}

int DataParser::get_string_property_array(
    const char * property_name, 
    std::string * string_array) 
{
    std::ifstream ifs;
    ifs.open(this->dat_filename.c_str());

    if (!ifs.good()) {

        fprintf(
            stderr, 
            "DataParser::get_string_property_array() : .scn file (%s) not found.\n", 
            this->dat_filename.c_str());

        ifs.close();

        return -1;
    }

    int n = 0;

    while (!ifs.eof()) {

        // read line of the .scn file
        char buf[MAX_CHARS_PER_LINE];
        ifs.getline(buf, MAX_CHARS_PER_LINE);

        // array to store memory addresses of the tokens in buf
        const char *token[MAX_TOKENS_PER_LINE] = {}; 

        // extract the first token, the property name
        token[0] = strtok(buf, PROPERTY_DELIMITER);

        if (token[0]) {

            if (strcmp(token[0], property_name) != 0) {
                // ...this isn't the property we're looking for. move along!
                continue;
            }

            for (n = 1; n < MAX_TOKENS_PER_LINE; ++n) {

                token[n] = strtok(0, VALUE_DELIMITER);

                if (!token[n]) {

                    break;

                } else {

                    string_array[n - 1] = std::string(token[n]);
                }
            }

            ifs.close();

            if (n > 0) {

                return 0;

            } else {

                fprintf(
                    stderr, 
                    "DataParser::get_string_property_array() : no tokens gathered\n");

                return -1;
            }
        }
    }

    ifs.close();

    fprintf(
        stderr, 
        "DataParser::get_string_property_array() : no property w/ name %s\n", 
        property_name);

    return -1;
}