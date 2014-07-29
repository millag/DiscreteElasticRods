#ifndef TRANSFORMTOOL_H
#define TRANSFORMTOOL_H

#include <ngl/Vec4.h>
#include <ngl/Mat4.h>

class TransformTool
{
public:
    TransformTool(ngl::Real phi = 0, ngl::Real theta = 0, ngl::Real scale = 1, ngl::Vec4 translation = ngl::Vec4(), ngl::Mat4 transform = ngl::Mat4()):
        m_phi(phi), m_theta(theta), m_scale(scale), m_translation(translation), m_transform(transform)
    {}

    void reset()
    {
        m_phi = m_theta = 0;
        m_scale = 1;
        m_translation.set(0,0,0,1);
        m_transform.identity();
    }

    ngl::Mat4 getTransform()
    {
        ngl::Mat4 rotX;
        ngl::Mat4 rotY;
        ngl::Mat4 scale;
        // create the rotation matrices
        rotX.rotateX(m_theta);
        rotY.rotateY(m_phi);
        scale.scale(m_scale, m_scale, m_scale);
        // multiply the rotations
        ngl::Mat4 m = rotY * rotX * scale;
        m_transform.translate(m_translation.m_x, m_translation.m_y, m_translation.m_z);
        return m * m_transform;
    }

    ngl::Real m_phi;
    ngl::Real m_theta;
    ngl::Real m_scale;
    ngl::Vec4 m_translation;
    ngl::Mat4 m_transform;
};

#endif // TRANSFORMTOOL_H
