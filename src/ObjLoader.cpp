#include "ObjLoader.h"
#include <iostream>
#include <fstream>


ObjLoader::ObjLoader(const char *filename)
{
    if (filename == NULL)
        return;

    loadFile(filename);
}

bool ObjLoader::loadFile(const char *filename) const
{
//     std::ifstream istream;
//     try { istream.open(filename);
//        if (istream.fail())
//        {
//            throw std::runtime_error("can't open file " + std::string(filename));
//        }

//        std::shared_ptr<Material> defaultMaterial(new Material);
//        curMaterial = defaultMaterial;
//        char line[MAX_LINE_LENGTH];
//        while (ifs.peek() != EOF)
//        {
//            ifs.getline(line, sizeof(line), '\n');
//            const char* token = line + strspn(line, " \t"); // ignore space and tabs
//            if (token[0] == 0)
//                continue;
//            if (token[0] == 'v' && isSep(token[1]))
//            {
//                v.push_back(getVec3f(token += 2));
//                continue;
//            }
//            if (!strncmp(token, "vn", 2) && isSep(token[2]))
//            {
//                vn.push_back(getVec3f(token += 3));
//                continue;
//            }
//            if (!strncmp(token, "vt", 2) && isSep(token[2]))
//            {
//                vt.push_back(getVec2f(token += 3));
//                continue;
//            }
//            if (token[0] == 'f' && isSep(token[1]))
//            {
//                parseSep(token += 1);
//                std::vector<Vertex> face;
//                while (token[0])
//                {
//                    face.push_back(getInt3(token));
//                    parseSepOpt(token);
//                }
//                curGroup.push_back(face);

//                continue;
//            }
//            /*! use material */
//            if (!strncmp(token, "usemtl", 6) && isSep(token[6]))
//            {
//                flushFaceGroup();
//                std::string name(parseSep(token += 6));
//                if (materials.find(name) == materials.end())
//                    curMaterial = defaultMaterial;
//                else
//                    curMaterial = materials[name];

//                continue;
//            }
//            /* load material library */
//            if (!strncmp(token, "mtllib", 6) && isSep(token[6]))
//            {
//                loadMTL(path + "/" + std::string(parseSep(token += 6)));
//                continue;
//            }
//        }
//     } catch (const std::exception &e)
//     {
//         std::cerr << e.what() << std::endl;
//     }
//     ifs.close();
}
