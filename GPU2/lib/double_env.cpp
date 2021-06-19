#include <stdlib.h>
#include <iostream>
#include <sstream>


extern "C" double double_env(const char* name) {
    double d = 0.0;
    const char* val = getenv(name);
    if (val) {
        std::istringstream ss(val);
        ss >> d;
        if (ss.fail()) {
            std::cerr << "WARNING: " << name << " \"" << val << "\" is invalid." << std::endl;
            d = 0.0;
        }
    }
    return d;
}
