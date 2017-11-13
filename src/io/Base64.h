#ifndef BASE64_H
#define BASE64_H

#include <string>

namespace IO {
void base64_encode(const char * , unsigned int len, std::string &out);
// void base64_decode(std::string const& s, std::string &out);
}

#endif

