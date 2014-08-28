#ifndef OBJPARSER_H
#define OBJPARSER_H

#include <vector>
#include "config.h"

class ObjLoader
{
public:

    ObjLoader(const char* filename = NULL);
    ~ObjLoader() { }

    bool loadFile(const char* filename);
    inline const std::vector<mg::Vec3D>& getVertices() const  { return m_v; }
    inline const std::vector<mg::Vec3D>& getNormals() const  { return m_vn; }
    inline const std::vector<mg::Vec2D>& getTexCoords() const  { return m_vt; }
    inline const std::vector<unsigned>& getVIndices() const  { return m_vindices; }
    inline const std::vector<unsigned>& getVNIndices() const  { return m_vnindices; }
    inline const std::vector<unsigned>& getVTIndices() const  { return m_vtindices; }

private:

    inline unsigned getVIdx(long int idx) { return (unsigned)(idx > 0 ? idx - 1 : (idx == 0 ? 0 : (long int)m_v.size() + idx)); }
    inline unsigned getVTIdx(long int idx) { return (unsigned)(idx > 0 ? idx - 1 : (idx == 0 ? 0 : (long int)m_vt.size() + idx)); }
    inline unsigned getVNIdx(long int idx) { return (unsigned)(idx > 0 ? idx - 1 : (idx == 0 ? 0 : (long int)m_vn.size() + idx)); }

    void parseVertex(std::string& s);
    void parseTexCoord(std::string& s);
    void parseNormal(std::string& s);
    void parseFace(std::string& s);

    std::vector<mg::Vec3D> m_v;
    std::vector<mg::Vec2D> m_vt;
    std::vector<mg::Vec3D> m_vn;
    std::vector<unsigned> m_vindices;
    std::vector<unsigned> m_vtindices;
    std::vector<unsigned> m_vnindices;
};

#endif // OBJPARSER_H
