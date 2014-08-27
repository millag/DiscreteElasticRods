#ifndef TRANSFORMTOOL_H
#define TRANSFORMTOOL_H

#include <ngl/Mat4.h>
#include "config.h"

/// TODO: need serious refactoring NO COMMENT :(
class TransformTool
{
public:
    TransformTool() { reset(); }

    void reset()
    {
        m_angle0 = m_angle1 = m_angle2 = 0;
        m_scale = 1;
        m_translation.set(0,0,0);
        m_matrix.identity();
    }

    ngl::Mat4 getMatrix()
    {
        ngl::Mat4 rotX;
        ngl::Mat4 rotY;
        ngl::Mat4 scale;
        // create the rotation matrices
        rotX.rotateX(m_angle1);
        rotY.rotateY(m_angle0);
        scale.scale(m_scale, m_scale, m_scale);
        // multiply the rotations
        ngl::Mat4 m = rotY * rotX * scale;
        m_matrix.translate(m_translation[0], m_translation[1], m_translation[2]);
        return m * m_matrix;
    }

    void setTransform(const mg::Matrix4D& m)
    {
        m_matrix.m_00 = m(0,0);
        m_matrix.m_01 = m(1,0);
        m_matrix.m_02 = m(2,0);
        m_matrix.m_03 = m(3,0);

        m_matrix.m_10 = m(0,1);
        m_matrix.m_11 = m(1,1);
        m_matrix.m_12 = m(2,1);
        m_matrix.m_13 = m(3,1);

        m_matrix.m_20 = m(0,2);
        m_matrix.m_21 = m(1,2);
        m_matrix.m_22 = m(2,2);
        m_matrix.m_23 = m(3,2);

        m_matrix.m_30 = m(0,3);
        m_matrix.m_31 = m(1,3);
        m_matrix.m_32 = m(2,3);
        m_matrix.m_33 = m(3,3);

        m_translation.set(m(0,3), m(1,3), m(2,3));
    }

    mg::Matrix4D getTransform()
    {
        ngl::Mat4 m = getMatrix();
        mg::Matrix4D transform;
        transform(0,0) = m.m_00;
        transform(1,0) = m.m_01;
        transform(2,0) = m.m_02;
        transform(3,0) = m.m_03;

        transform(0,1) = m.m_10;
        transform(1,1) = m.m_11;
        transform(2,1) = m.m_12;
        transform(3,1) = m.m_13;

        transform(0,2) = m.m_20;
        transform(1,2) = m.m_21;
        transform(2,2) = m.m_22;
        transform(3,2) = m.m_23;

        transform(0,3) = m.m_30;
        transform(1,3) = m.m_31;
        transform(2,3) = m.m_32;
        transform(3,3) = m.m_33;

        return transform;
    }

    mg::Real m_angle0;
    mg::Real m_angle1;
    mg::Real m_angle2;
    mg::Real m_scale;
    mg::Vec3D m_translation;
    ngl::Mat4 m_matrix;
};

class Transformation
{
public:
    Transformation() { reset(); }

    void reset() { m_matrix.identity(); }

    const ngl::Mat4& getMatrix() const  { return m_matrix; }
    void setMatrix(const mg::Matrix4D& m)
    {
        m_matrix.m_00 = m(0,0);
        m_matrix.m_01 = m(1,0);
        m_matrix.m_02 = m(2,0);
        m_matrix.m_03 = m(3,0);

        m_matrix.m_10 = m(0,1);
        m_matrix.m_11 = m(1,1);
        m_matrix.m_12 = m(2,1);
        m_matrix.m_13 = m(3,1);

        m_matrix.m_20 = m(0,2);
        m_matrix.m_21 = m(1,2);
        m_matrix.m_22 = m(2,2);
        m_matrix.m_23 = m(3,2);

        m_matrix.m_30 = m(0,3);
        m_matrix.m_31 = m(1,3);
        m_matrix.m_32 = m(2,3);
        m_matrix.m_33 = m(3,3);
    }

private:
    ngl::Mat4 m_matrix;
};


#endif // TRANSFORMTOOL_H
