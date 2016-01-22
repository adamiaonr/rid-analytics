#include "dataparser.h"

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

        fprintf(stderr, "DataParser::get_int_property_array() : got property = %s\n", token[0]);

        if (token[0]) {

            for (n = 0; n < MAX_TOKENS_PER_LINE; ++n) {

                token[n] = strtok(0, VALUE_DELIMITER);

                if (!token[n]) {

                    break;

                } else {

                    fprintf(stderr, "DataParser::get_int_property_array() : got value = %s\n", token[n]);

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

int DataParser::get_double_property_value(const char * property_name, double &double_value) 
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
    double * double_array) 
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

                    double_array[n] = atof(token[n]);
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