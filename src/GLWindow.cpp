#include "GLWindow.h"
#include <iostream>
#include <ngl/Camera.h>
#include <ngl/Colour.h>
#include <ngl/Light.h>
#include <ngl/Mat4.h>
#include <ngl/Transformation.h>
#include <ngl/TransformStack.h>
#include <ngl/Material.h>
#include <ngl/NGLInit.h>
#include <ngl/Obj.h>
#include <ngl/Random.h>
#include <ngl/VAOPrimitives.h>
#include <ngl/ShaderLib.h>
#include "Utils.h"

//----------------------------------------------------------------------------------------------------------------------
GLWindow::GLWindow(const QGLFormat _format, QWidget *_parent ) : QGLWidget( _format, _parent )
{
    m_rotate = false;
    m_zoom = false;
    m_selection = false;

    m_scene = NULL;
    m_strandVAO = NULL;

    // re-size the widget to that of the parent (in this case the GLFrame passed in on construction)
    this->resize(_parent->size());
    // set this widget to have the initial keyboard focus
    setFocus();
}

GLWindow::~GLWindow()
{
    std::cout<<"Shutting down NGL, removing VAO's and Shaders\n";

    delete m_cam;

    ngl::NGLInit *init = ngl::NGLInit::instance();
    init->NGLQuit();
}

// This virtual function is called once before the first call to paintGL() or resizeGL(),
//and then once whenever the widget has been assigned a new QGLContext.
// This function should set up any required OpenGL context rendering flags, defining display lists, etc.

//----------------------------------------------------------------------------------------------------------------------
void GLWindow::initializeGL()
{
    m_cameraTransform.m_angle0 = 0;
    m_cameraTransform.m_angle1 = 10;
    m_cameraTransform.m_translation.set(0,0,-17);

    ngl::NGLInit::instance();
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);			   // Grey Background
    // enable depth testing for drawing
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LINE_SMOOTH);
    // Now we will create a basic Camera from the graphics library
    // This is a static camera so it only needs to be set once
    // First create Values for the camera position
    ngl::Vec3 from(0,0,1);
    ngl::Vec3 to(0,0,0);
    ngl::Vec3 up(0,1,0);
    ngl::NGLInit::instance();
    m_cam = new ngl::Camera(from,to,up);
    // set the shape using FOV 45 Aspect Ratio based on Width and Height
    // The final two are near and far clipping planes of 0.5 and 10
    m_cam->setShape(45,(float)720.0/576.0,0.5,150);
    // now to load the shader and set the values
    // grab an instance of shader manager
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();

    // we are creating a shader called Phong
    shader->createShaderProgram("Phong");
    // now we are going to create empty shaders for Frag and Vert
    shader->attachShader("PhongVertex",ngl::VERTEX);
    shader->attachShader("PhongFragment",ngl::FRAGMENT);
    // attach the source
    shader->loadShaderSource("PhongVertex","shaders/PhongVertex.glsl");
    shader->loadShaderSource("PhongFragment","shaders/PhongFragment.glsl");
    // compile the shaders
    shader->compileShader("PhongVertex");
    shader->compileShader("PhongFragment");
    // add them to the program
    shader->attachShaderToProgram("Phong","PhongVertex");
    shader->attachShaderToProgram("Phong","PhongFragment");
    // now bind the shader attributes for most NGL primitives we use the following
    // layout attribute 0 is the vertex data (x,y,z)
    shader->bindAttribute("Phong",0,"inVert");
    // attribute 1 is the UV data u,v (if present)
    shader->bindAttribute("Phong",1,"inUV");
    // attribute 2 are the normals x,y,z
    shader->bindAttribute("Phong",2,"inNormal");

    // now we have associated this data we can link the shader
    shader->linkProgramObject("Phong");
    // and make it active ready to load values
    (*shader)["Phong"]->use();
    shader->setShaderParam1i("Normalize",1);

    // now pass the modelView and projection values to the shader
    // the shader will use the currently active material and light0 so set them
    ngl::Material m(ngl::SILVER);
    m.loadToShader("material");
    ngl::Light light(ngl::Vec3(0,0,0),ngl::Colour(1,1,1,1),ngl::Colour(1,1,1,1),ngl::POINTLIGHT);
    // now create our light this is done after the camera so we can pass the
    // transpose of the projection matrix to the light to do correct eye space
    // transformations
    ngl::Mat4 iv=m_cam->getViewMatrix();
    iv.transpose();
    light.setTransform(iv);
    light.setAttenuation(1,0,0);
    light.enable();
    // load these values to the shader as well
    light.loadToShader("light");


    shader->createShaderProgram("Tube");

    shader->attachShader("TubeVertex", ngl::VERTEX);
    shader->attachShader("TubeTCS", ngl::TESSCONTROL);
    shader->attachShader("TubeTES", ngl::TESSEVAL);
    shader->attachShader("TubeFragment", ngl::FRAGMENT);
    shader->loadShaderSource("TubeVertex", "shaders/TubeVert.glsl");
    shader->loadShaderSource("TubeTCS", "shaders/TubeTCS.glsl");
    shader->loadShaderSource("TubeTES", "shaders/TubeTES.glsl");
    shader->loadShaderSource("TubeFragment", "shaders/TubeFrag.glsl");

    shader->compileShader("TubeVertex");
    shader->compileShader("TubeTCS");
    shader->compileShader("TubeTES");
    shader->compileShader("TubeFragment");
    shader->attachShaderToProgram("Tube","TubeVertex");
    shader->attachShaderToProgram("Tube","TubeTCS");
    shader->attachShaderToProgram("Tube","TubeTES");
    shader->attachShaderToProgram("Tube","TubeFragment");

    shader->bindAttribute("Tube",0,"inVert");
    shader->bindAttribute("Tube",1,"inNormal");

    shader->linkProgramObject("Tube");
    (*shader)["Tube"]->use();
    light.loadToShader("light");
    m.loadToShader("material");


    shader->createShaderProgram("DebugRod");

    shader->attachShader("DebugGeometry", ngl::GEOMETRY);
    shader->attachShader("DebugVertex", ngl::VERTEX);
    shader->attachShader("DebugFragment", ngl::FRAGMENT);
    shader->loadShaderSource("DebugGeometry", "shaders/DebugGeom.glsl");
    shader->loadShaderSource("DebugVertex", "shaders/DebugVert.glsl");
    shader->loadShaderSource("DebugFragment", "shaders/DebugFrag.glsl");

    shader->compileShader("DebugGeometry");
    shader->compileShader("DebugVertex");
    shader->compileShader("DebugFragment");
    shader->attachShaderToProgram("DebugRod","DebugGeometry");
    shader->attachShaderToProgram("DebugRod","DebugVertex");
    shader->attachShaderToProgram("DebugRod","DebugFragment");

    shader->bindAttribute("DebugRod",0,"inVert");
    shader->bindAttribute("DebugRod",1,"inKB");
    shader->bindAttribute("DebugRod",2,"inM1");
    shader->bindAttribute("DebugRod",3,"inM2");
    shader->bindAttribute("DebugRod",4,"inForce");

    shader->linkProgramObject("DebugRod");
    (*shader)["DebugRod"]->use();
    light.loadToShader("light");
    m.loadToShader("material");


    ngl::VAOPrimitives *prim=ngl::VAOPrimitives::instance();
    prim->createLineGrid("grid", 10, 10, 10);

    if (m_scene == NULL)
    {
        return;
    }

//    build VAO for each mesh in scene and store references in m_VAOList
    buildVAOs(m_scene->getMeshes(), m_VAOList);
}

//----------------------------------------------------------------------------------------------------------------------
//This virtual function is called whenever the widget has been resized.
// The new size is passed in width and height.
void GLWindow::resizeGL( int _w, int _h)
{
    glViewport(0, 0, _w, _h);
    m_cam->setShape(45, (float)_w/_h, 0.5, 150);
}

void GLWindow::loadMatricesToShader()
{
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();

    ngl::Mat4 MV;
    ngl::Mat4 MVP;
    ngl::Mat3 normalMatrix;
    ngl::Mat4 M;
    M = m_transform.getMatrix();
    MV = M * m_cameraTransform.getMatrix() * m_cam->getViewMatrix();
    MVP = MV * m_cam->getProjectionMatrix();
    normalMatrix = MV;
    normalMatrix.inverse();
    shader->setShaderParamFromMat4("MV",MV);
    shader->setShaderParamFromMat4("MVP",MVP);
    shader->setShaderParamFromMat3("normalMatrix",normalMatrix);
    shader->setShaderParamFromMat4("M",M);
}

void GLWindow::loadMatricesToHairShader()
{
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();

    ngl::Mat4 M;
    ngl::Mat4 MV;
    ngl::Mat4 MVP;
    M = m_transform.getMatrix();
    MV = M  * m_cameraTransform.getMatrix() * m_cam->getViewMatrix();
    MVP = MV * m_cam->getProjectionMatrix();
    shader->setShaderParamFromMat4("MVP",MVP);
    shader->setShaderParamFromMat4("M",M);
}

void GLWindow::loadMatricesToHairShader2()
{
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();

    ngl::Mat4 MV;
    ngl::Mat4 MVP;
    ngl::Mat3 normalMatrix;
    ngl::Mat4 M;
    M = m_transform.getMatrix();
    MV = M * m_cameraTransform.getMatrix() * m_cam->getViewMatrix();
    MVP = MV * m_cam->getProjectionMatrix();
    normalMatrix = MV;
    normalMatrix.inverse();
    shader->setShaderParamFromMat4("MV",MV);
    shader->setShaderParamFromMat4("MVP",MVP);
    shader->setShaderParamFromMat3("normalMatrix",normalMatrix);
}

//----------------------------------------------------------------------------------------------------------------------
//This virtual function is called whenever the widget needs to be painted.
// this is our main drawing routine
void GLWindow::paintGL()
{
//    clear the screen and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

//    grab an instance of the shader manager
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();
//    get an instance of the VBO primitives for drawing
    ngl::VAOPrimitives *prim=ngl::VAOPrimitives::instance();

//    draw grid
//    shader->use("Phong");
//    m_transform.reset();
//    loadMatricesToShader();
//    prim->draw("grid");

    if (m_scene == NULL)
    {
        return;
    }

    ngl::VertexArrayObject* vao;
    const std::vector<RenderObject*>& roList = m_scene->getRenderObjects();
    typedef std::vector<RenderObject*>::const_iterator ROIter;
    for (ROIter it = roList.begin(); it != roList.end(); ++it)
    {
        RenderObject* ro = (*it);
        m_transform.reset();
        {
//            place object in world
            shader->use("Phong");
            m_transform.setMatrix(ro->getTransform());
            loadMatricesToShader();

//            actual draw
            vao =  m_VAOList.at(ro->getMeshId());
            if (vao == NULL)
                continue;

            vao->bind();
            vao->draw();
            vao->unbind();
        }
    }

    if (m_scene->getHairById(0) == NULL)
    {
        return;
    }

#ifndef DBUGG
    shader->use("Tube");
    shader->setShaderParam4f("Col1", 0.47, 0.4, 0.0, 1.0);
    shader->setShaderParam4f("Col2", 0.0, 0.0, 0.0, 1.0);
    shader->setShaderParam1f("Radius", 0.1);
    glPatchParameteri(GL_PATCH_VERTICES, 4);
    // load transform stack
    m_transform.reset();
    loadMatricesToHairShader2();
#endif

#ifdef DBUGG
    glLineWidth(0.05);
    shader->use("DebugRod");
    shader->setShaderParam4f("Colour", 0.8, 0.8, 0.0, 1.0);
    // load transform stack
    m_transform.reset();
    loadMatricesToHairShader();
#endif

    const std::vector<ElasticRod*>& strands = m_scene->getHairById(0)->m_strands;
//    const std::vector<ElasticRod*>& strands = m_scene->getStrands();
    typedef std::vector<ElasticRod*>::const_iterator SIter;
    for (SIter it = strands.begin(); it != strands.end(); ++it)
    {
        drawHairStrand(**it);
    }
}





//----------------------------------------------------------------------------------------------------------------------
void GLWindow::mouseMoveEvent ( QMouseEvent * _event )
{
    TransformTool& ttool = (m_selection)? m_selectionTransform : m_cameraTransform;

    // note the method buttons() is the button state when event was called
    // this is different from button() which is used to check which button was
    // pressed when the mousePress/Release event is generated
    if(m_rotate && _event->buttons() == Qt::LeftButton)
    {
        mg::Real dx = _event->x() - m_mouseX;
        mg::Real dy = _event->y() - m_mouseY;
        ttool.m_angle0 = fmod(ttool.m_angle0 + dx, 360.0);
        ttool.m_angle1 = fmod(ttool.m_angle1 + dy, 360.0);

        m_mouseX = _event->x();
        m_mouseY = _event->y();
        // re-draw GL
        update();
    }

    if(m_zoom && _event->buttons() == Qt::RightButton)
    {
        mg::Real dx = 2 * (mg::Real)(_event->x() - m_mouseX) / width();
        mg::Real dy = 2 * (mg::Real)(_event->y() - m_mouseY) / height();
        ngl::Vec4 lookDir = (m_cam->getLook() - m_cam->getEye());
        mg::Vec3D dirz(lookDir.m_x, lookDir.m_y, lookDir.m_z);

        if (m_selection)
        {
            dirz = mg::transform_vector(m_cameraTransform.getTransform().inverse(), dirz);
            dirz.normalize();
            mg::Vec3D dirx = mg::cross(mg::EY, dirz);
            dirx.normalize();
            mg::Vec3D diry = mg::cross(dirz, dirx);
            diry.normalize();
            ttool.m_translation += dirx * dx;
            ttool.m_translation += diry * dy;
        } else
        {
            ttool.m_translation = ttool.m_translation + dirz * dx;
        }

        m_mouseX = _event->x();
        m_mouseY = _event->y();
        // re-draw GL
        update();
    }
}


//----------------------------------------------------------------------------------------------------------------------
void GLWindow::mousePressEvent ( QMouseEvent * _event  )
{
    // this method is called when the mouse button is pressed in this case we
    // store the value where the maouse was clicked (x,y) and set the Rotate flag to true

    m_mouseX = _event->x();
    m_mouseY = _event->y();

    switch (_event->button()) {
        case Qt::LeftButton:
            m_rotate = true;
            break;
        case Qt::RightButton:
            m_zoom = true;
            break;

        default:
            break;
    }
}

//----------------------------------------------------------------------------------------------------------------------
void GLWindow::mouseReleaseEvent ( QMouseEvent * _event )
{
    // this event is called when the mouse button is released
    // we then set Rotate to false

    switch (_event->button()) {
        case Qt::LeftButton:
            m_rotate = false;
            break;
        case Qt::RightButton:
            m_zoom = false;
            break;

        default:
            break;
    }
}

// ============================ utility functions ==========================

void GLWindow::buildVAOs(const std::vector<Mesh*>& meshList, std::vector<ngl::VertexArrayObject*>& o_VAOList) const
{
    assert(o_VAOList.size() == 0);
    o_VAOList.resize(meshList.size(), NULL);

    unsigned i = 0;
    for (std::vector<Mesh*>::const_iterator it = meshList.begin(); it != meshList.end(); ++it)
    {
        if ((*it)->getNPrimitives() == 0)
        {
            ++i;
            continue;
        }

        // create a vao as a series of GL_TRIANGLES
        ngl::VertexArrayObject* vao = ngl::VertexArrayObject::createVOA(GL_TRIANGLES);
        o_VAOList[i++] = vao;
        // feed data to GPU
        feedVAO(*(*it), *vao);
    }
}

void GLWindow::feedVAO(const Mesh& mesh, ngl::VertexArrayObject& o_vao) const
{
    o_vao.bind();
    // in this case we are going to set our data as the vertices above
    // now we set the attribute pointer to be 0 (as this matches vertIn in our shader)
    o_vao.setIndexedData(mesh.m_vertices.size() * sizeof(mg::Vec3D),
                         mesh.m_vertices[0][0],
                         mesh.m_vindices.size(),
                         &mesh.m_vindices[0],
                         GL_UNSIGNED_INT);
    o_vao.setVertexAttributePointer(0, 3, GL_FLOAT, 0, 0);

    // now we set the attribute pointer to be 2 (as this matches normalIn in our shader)
    o_vao.setIndexedData(mesh.m_normals.size() * sizeof(mg::Vec3D),
                         mesh.m_normals[0][0],
                         mesh.m_vindices.size(),
                         &mesh.m_vindices[0],
                         GL_UNSIGNED_INT);
    o_vao.setVertexAttributePointer(2, 3, GL_FLOAT, 0, 0);

    o_vao.setNumIndices(mesh.m_vindices.size());

    // now unbind
    o_vao.unbind();
}


#ifdef DBUGG
void GLWindow::drawHairStrand(const ElasticRod& strand)
{
    if (m_strandVAO != NULL)
    {
        m_strandVAO->bind();

        m_strandVAO->updateIndexedData(0, strand.m_ppos.size() * sizeof(mg::Vec3D),
                                       strand.m_ppos[0][0]);
        m_strandVAO->setVertexAttributePointer(0, 3, GL_FLOAT, 0, 0);

        m_strandVAO->updateIndexedData(1, strand.m_kb.size() * sizeof(mg::Vec3D),
                                       strand.m_kb[0][0]);
        m_strandVAO->setVertexAttributePointer(1, 3, GL_FLOAT, 0, 0);

        m_strandVAO->updateIndexedData(2, strand.m_m1.size() * sizeof(mg::Vec3D),
                                       strand.m_m1[0][0]);
        m_strandVAO->setVertexAttributePointer(2, 3, GL_FLOAT, 0, 0);

        m_strandVAO->updateIndexedData(3, strand.m_m2.size() * sizeof(mg::Vec3D),
                                       strand.m_m2[0][0]);
        m_strandVAO->setVertexAttributePointer(3, 3, GL_FLOAT, 0, 0);

        m_strandVAO->updateIndexedData(4, strand.m_elasticForce.size() * sizeof(mg::Vec3D),
                                       strand.m_elasticForce[0][0]);
        m_strandVAO->setVertexAttributePointer(4, 3, GL_FLOAT, 0, 0);

        m_strandVAO->draw();
        m_strandVAO->unbind();
        return;
    }

    std::vector<unsigned> indices(strand.m_ppos.size() + 2);
    unsigned i = 0;
    indices[i] = 0;
    for (i = 0; i < strand.m_ppos.size(); ++i)
    {
        indices[i + 1] = i;
    }
    indices[i + 1] = i - 1;

    m_strandVAO = ngl::VertexArrayObject::createVOA(GL_LINE_STRIP_ADJACENCY);
    m_strandVAO->bind();

    m_strandVAO->setIndexedData(strand.m_ppos.size() * sizeof(mg::Vec3D),
                                strand.m_ppos[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(0, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setIndexedData(strand.m_kb.size() * sizeof(mg::Vec3D),
                                strand.m_kb[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(1, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setIndexedData(strand.m_m1.size() * sizeof(mg::Vec3D),
                                strand.m_m1[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(2, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setIndexedData(strand.m_m2.size() * sizeof(mg::Vec3D),
                                strand.m_m2[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(3, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setIndexedData(strand.m_elasticForce.size() * sizeof(mg::Vec3D),
                                strand.m_elasticForce[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(4, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setNumIndices(indices.size());

    m_strandVAO->draw();
    m_strandVAO->unbind();
}
#endif

#ifndef DBUGG
void GLWindow::drawHairStrand(const ElasticRod& strand)
{
    if (m_strandVAO != NULL)
    {
        m_strandVAO->bind();

        m_strandVAO->updateIndexedData(0, strand.m_ppos.size() * sizeof(mg::Vec3D),
                                       strand.m_ppos[0][0]);
        m_strandVAO->setVertexAttributePointer(0, 3, GL_FLOAT, 0, 0);

        m_strandVAO->updateIndexedData(1, strand.m_m1.size() * sizeof(mg::Vec3D),
                                       strand.m_m1[0][0]);
        m_strandVAO->setVertexAttributePointer(1, 3, GL_FLOAT, 0, 0);

        m_strandVAO->draw();
        m_strandVAO->unbind();
        return;
    }

    std::vector<unsigned> indices( (strand.m_ppos.size() - 1) * 4 );
    int nPoints = (int)strand.m_ppos.size();
    for (int i = -1; i < nPoints - 2; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            assert( (i + 1) >= 0 && (4*(i + 1) + j) < (int)indices.size());
            indices[4*(i + 1) + j] = std::min(std::max(0, i + j), nPoints - 1);
        }
    }

    m_strandVAO = ngl::VertexArrayObject::createVOA(GL_PATCHES);
    m_strandVAO->bind();

    m_strandVAO->setIndexedData(strand.m_ppos.size() * sizeof(mg::Vec3D),
                                strand.m_ppos[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(0, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setIndexedData(strand.m_m1.size() * sizeof(mg::Vec3D),
                                strand.m_m1[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(1, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setNumIndices(indices.size());

    m_strandVAO->draw();
    m_strandVAO->unbind();
}
#endif
