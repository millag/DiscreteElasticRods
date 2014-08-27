#ifndef OBJPARSER_H
#define OBJPARSER_H

#include <vector>
#include "config.h"

class ObjLoader
{
public:

    ObjLoader(const char* filename = NULL);
    ~ObjLoader() { }

    bool loadFile(const char* filename) const;
    const std::vector<mg::Vec3D>& getVertices() const  { return m_v; }
    const std::vector<mg::Vec3D>& getNormals() const  { return m_vn; }
    const std::vector<mg::Vec2D>& getTexCoords() const  { return m_vt; }
    const std::vector<unsigned>& getVIndices() const  { return m_vindices; }
    const std::vector<unsigned>& getVNIndices() const  { return m_vnindices; }
    const std::vector<unsigned>& getVTIndices() const  { return m_vtindices; }

private:
//    Vertex getInt3(const char*& token);
    inline int fix_v(int idx) { return(idx > 0 ? idx - 1 : (idx == 0 ? 0 : (int)m_v.size() + idx)); }
    inline int fix_vt(int idx) { return(idx > 0 ? idx - 1 : (idx == 0 ? 0 : (int)m_vt.size() + idx)); }
    inline int fix_vn(int idx) { return(idx > 0 ? idx - 1 : (idx == 0 ? 0 : (int)m_vn.size() + idx)); }
    std::vector<mg::Vec3D> m_v;
    std::vector<mg::Vec3D> m_vn;
    std::vector<mg::Vec2D> m_vt;
    std::vector<unsigned> m_vindices;
    std::vector<unsigned> m_vnindices;
    std::vector<unsigned> m_vtindices;

//    void flushFaceGroup();
//    uint32_t getVertex(std::map<Vertex, uint32_t>&, std::vector<Vec3f>&, std::vector<Vec3f>&, std::vector<Vec2f>&, const Vertex&);
//    std::vector<std::shared_ptr<Primitive> > model;
};

#endif // OBJPARSER_H
