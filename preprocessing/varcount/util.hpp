#ifndef UTIL_HPP
#define UTIL_HPP

#include <vector>
#include <cstdio>
#include <cstring>
#include <string>


std::vector<std::string> parse_ids(const char* s) {
    std::vector<std::string> ss;
    char* s__ = (char*) malloc(sizeof(char) * (strlen(s) + 1));
    strcpy(s__, s);
    char* s_ = std::strtok(s__, ";");
    while (s_ != NULL) {
        ss.push_back(s_);
        s_ = std::strtok(NULL, ";");
    }
    free(s__);
    return ss;
}

#endif
