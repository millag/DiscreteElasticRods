#include "ObjLoader.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include "BasicParser.h"

ObjLoader::ObjLoader(const char *filename)
{
    loadFile(filename);
}

bool ObjLoader::loadFile(const char *filename)
{
    if (filename == NULL)
        return false;

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
        BasicParser::trim(line);
        token = BasicParser::parseToken(line);
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

    ifs.close();
    return true;
}

void ObjLoader::parseVertex(std::string& s)
{
    mg::Vec3D v;
    BasicParser::parseVec3D(s, v);
    m_v.push_back(v);
}

void ObjLoader::parseTexCoord(std::string& s)
{
    mg::Vec2D vt;
    BasicParser::parseVec2D(s, vt);
    m_vt.push_back(vt);
}

void ObjLoader::parseNormal(std::string& s)
{
    mg::Vec3D vn;
    BasicParser::parseVec3D(s, vn);
    m_vn.push_back(vn);
}

void ObjLoader::parseFace(std::string& s)
{
    BasicParser::trim(s);
    std::vector<std::string> tokens;
    BasicParser::split(s, ' ', tokens);
    for (unsigned i = 0; i < tokens.size(); ++i)
    {
        std::vector<std::string> indices;
        BasicParser::split(tokens[i], '/', indices);
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
