#include "ObjLoader.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <cstdlib>
#include <cassert>

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

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
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

static inline float parseFloat(const std::string &s)
{
    return std::strtod(s.c_str(), NULL);
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

ObjLoader::ObjLoader(const char *filename)
{
    if (filename == NULL)
        return;

    loadFile(filename);
}

bool ObjLoader::loadFile(const char *filename)
{
    std::ifstream ifs;
    ifs.open(filename, std::ios::in);
    if (ifs.fail())
    {
        std::cerr << "can't open file " + std::string(filename) << std::endl;
        ifs.close();
        return false;
    }

    std::string line;
    std::string token;
    while (ifs.peek() != EOF)
    {
        std::getline(ifs, line, '\n');
        trim(line);
        token = parseToken(line);
        if (token.empty())
        {
            continue;
        }

        if (token.compare("v") == 0)
        {
            line.erase(0, token.length());
            parseVertex(line);
            continue;
        }
        if (token.compare("vt") == 0)
        {
            line.erase(0, token.length());
            parseTexCoord(line);
            continue;
        }
        if (token.compare("vn") == 0)
        {
            line.erase(0, token.length());
            parseNormal(line);
            continue;
        }
        if (token.compare("f") == 0)
        {
            line.erase(0, token.length());
            parseFace(line);
            continue;
        }
    }
    return true;
}

void ObjLoader::parseVertex(std::string& s)
{
    mg::Vec3D v;
    parseVec3D(s, v);
    m_v.push_back(v);
}

void ObjLoader::parseTexCoord(std::string& s)
{
    mg::Vec2D vt;
    parseVec2D(s, vt);
    m_vt.push_back(vt);
}

void ObjLoader::parseNormal(std::string& s)
{
    mg::Vec3D vn;
    parseVec3D(s, vn);
    m_vn.push_back(vn);
}

void ObjLoader::parseFace(std::string& s)
{
    trim(s);
    std::vector<std::string> tokens;
    split(s, ' ', tokens);
    for (unsigned i = 0; i < tokens.size(); ++i)
    {
        std::vector<std::string> indices;
        split(tokens[i], '/', indices);
        assert(indices.size() > 0);

        m_vindices.push_back( getVIdx(std::strtol(indices[0].c_str(), NULL, 10)) );
        if (indices.size() > 1 && !indices[1].empty())
        {
            m_vtindices.push_back( getVIdx(std::strtol(indices[1].c_str(), NULL, 10)) );
        }
        if (indices.size() > 2)
        {
            m_vnindices.push_back( getVIdx(std::strtol(indices[2].c_str(), NULL, 10)) );
        }
    }
}
