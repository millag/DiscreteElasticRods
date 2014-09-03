#ifndef BASICPARSER_H
#define BASICPARSER_H

#include "config.h"

#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <cstdlib>

class BasicParser
{
public:

    // Citation Begin
    // code borrowed from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring

    // trim from start
    static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }

    // trim from end
    static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    // trim from both ends
    static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
    }

    // Citation End

    // Citation Begin
    // code borrowed from http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c

    static inline std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
        return elems;
    }

    // Citation End


    static inline std::string parseToken(const std::string &s)
    {
        std::string::const_iterator it = std::find_if(s.begin(), s.end(), std::ptr_fun<int, int>(std::isspace));
        return s.substr(0, it - s.begin());
    }

    static inline std::string parseWord(const std::string &s)
    {
        std::string::const_iterator it = std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isalpha)));
        return s.substr(0, it - s.begin());
    }

    static inline double parseDouble(const std::string &s)
    {
        return std::strtod(s.c_str(), NULL);
    }

    static inline long int parseInt(const std::string &s)
    {
        return std::strtol(s.c_str(), NULL, 10);
    }

    static inline void parseVec3D(std::string& s, mg::Vec3D &o_v)
    {
        char* end;
        o_v[0] = std::strtod(s.c_str(), &end);
        o_v[1] = std::strtod(end, &end);
        o_v[2] = std::strtod(end, &end);
    }

    static inline void parseVec2D(std::string& s, mg::Vec2D &o_v)
    {
        char* end;
        o_v[0] = std::strtod(s.c_str(), &end);
        o_v[1] = std::strtod(end, &end);
    }

    static inline void parseMatrix4D(std::string& s, mg::Matrix4D &o_m)
    {
        char* end;
        mg::Real val_x = std::strtod(s.c_str(), &end);
        mg::Real val_y = std::strtod(end, &end);
        mg::Real val_z = std::strtod(end, &end);
        mg::matrix_rotation_euler(o_m, mg::rad(val_x), mg::rad(val_y), mg::rad(val_z), mg::euler_order_xyz);

        val_x = std::strtod(end, &end);
        val_y = std::strtod(end, &end);
        val_z = std::strtod(end, &end);
        mg::Matrix4D scale;
        mg::matrix_scale(scale, val_x, val_y, val_z);
        o_m *= scale;

        val_x = std::strtod(end, &end);
        val_y = std::strtod(end, &end);
        val_z = std::strtod(end, &end);
        mg::matrix_set_translation(o_m, val_x, val_y, val_z);
    }

    static inline void parseIntList(std::string& s, std::vector<unsigned> &o_list)
    {
        const char* start = s.c_str();
        char* end;
        for (unsigned i = 0; i < o_list.size(); ++i)
        {
            o_list[i] = std::strtod(start, &end);
            start = end;
        }
    }

};

#endif // BASICPARSER_H
