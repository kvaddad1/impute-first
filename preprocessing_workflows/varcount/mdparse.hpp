#ifndef PARSE_MD_CPP
#define PARSE_MD_CPP

#include <vector>
#include <string>
#include <cstdio>

#define MD_END -1
#define MD_SNP 0
#define MD_DEL 1

struct MDPos {
    uint32_t p = 0;
    std::string str = ""; // maybe we can get away with char* instead?
    int st = MD_END;
};

// x
std::vector<MDPos>& md_parse(const char* md, std::vector<MDPos>& mdv) {
    mdv.clear();
    const char* s = &md[0];
    uint32_t pos = 0;
    while (*s) {
        MDPos mdp;
        while (*s && std::isdigit(s[0])) {
            mdp.p = (mdp.p*10) + (s[0] - '0');
            ++s;
        } 
        mdp.p += pos;
        pos = mdp.p;
        while (*s && !std::isdigit(s[0])) {
            if (mdp.st != MD_DEL) {
                mdp.st = MD_SNP;
            }
            if (s[0] == '^') mdp.st = MD_DEL;
            else mdp.str += s[0]; 
            ++s;
        }
        pos += mdp.str.size();
        mdv.push_back(mdp);
    }
    return mdv;
}

std::vector<MDPos> md_parse(const char* md) {
    std::vector<MDPos> v;
    return md_parse(md, v);
}

#endif

