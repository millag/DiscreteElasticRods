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
#include <boost/foreach.hpp>

#include "Utils.h"

void buildVAOs(const std::vector<Mesh*>& meshList, std::vector<ngl::VertexArrayObject*>& o_VAOList);
void feedVAO(const Mesh &mesh, ngl::VertexArrayObject& o_vao);

//----------------------------------------------------------------------------------------------------------------------
GLWindow::GLWindow(const QGLFormat _format, int _timer, QWidget *_parent ) : QGLWidget( _format, _parent )
{
    // re-size the widget to that of the parent (in this case the GLFrame passed in on construction)
    this->resize(_parent->size());
    // set this widget to have the initial keyboard focus
    setFocus();

    // Now set the initial GLWindow attributes to default values
    // Roate is false
    m_selectedObject = NULL;
    m_rotate = false;
    m_zoom = false;
    m_pan = false;

    m_strandVAO = NULL;

    m_timerValue = _timer;
    startSimTimer();
}

GLWindow::~GLWindow()
{
    std::cout<<"Shutting down NGL, removing VAO's and Shaders\n";
    ngl::NGLInit *init = ngl::NGLInit::instance();

    if (m_scene != NULL)
    {
        delete m_cam;
    }

    init->NGLQuit();
}

// This virtual function is called once before the first call to paintGL() or resizeGL(),
//and then once whenever the widget has been assigned a new QGLContext.
// This function should set up any required OpenGL context rendering flags, defining display lists, etc.

//----------------------------------------------------------------------------------------------------------------------
void GLWindow::initializeGL()
{
    m_cameraTransform.m_phi = 0;
    m_cameraTransform.m_theta = 30;
    m_cameraTransform.m_translation.set(0,0,-6);

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
    ngl::Light light(ngl::Vec3(3,3,0),ngl::Colour(1,1,1,1),ngl::Colour(1,1,1,1),ngl::POINTLIGHT);
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


    shader->createShaderProgram("Colour");

    shader->attachShader("ColourGeometry", ngl::GEOMETRY);
    shader->attachShader("ColourVertex", ngl::VERTEX);
    shader->attachShader("ColourFragment", ngl::FRAGMENT);
    shader->loadShaderSource("ColourGeometry", "shaders/ColourGeom.glsl");
    shader->loadShaderSource("ColourVertex", "shaders/ColourVert.glsl");
    shader->loadShaderSource("ColourFragment", "shaders/ColourFrag.glsl");

    shader->compileShader("ColourGeometry");
    shader->compileShader("ColourVertex");
    shader->compileShader("ColourFragment");
    shader->attachShaderToProgram("Colour","ColourGeometry");
    shader->attachShaderToProgram("Colour","ColourVertex");
    shader->attachShaderToProgram("Colour","ColourFragment");

    shader->bindAttribute("Colour",0,"inVert");

    shader->linkProgramObject("Colour");
    (*shader)["Colour"]->use();
    shader->setShaderParam1i("Normalize",1);
    light.loadToShader("light");


    ngl::VAOPrimitives *prim=ngl::VAOPrimitives::instance();
    prim->createSphere("sphere",1.0,20);

    prim->createLineGrid("grid", 10, 10, 10);

    // build VAO for each mesh in scene and store references in m_VAOList
    if (m_scene)
    {
        buildVAOs(m_scene->getMeshes(), m_VAOList);
    }
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
    MV = M * m_cameraTransform.getTransform() * m_cam->getViewMatrix();
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
    MV = M * m_cameraTransform.getTransform() * m_cam->getViewMatrix();
    MVP = MV * m_cam->getProjectionMatrix();
    shader->setShaderParamFromMat4("MVP",MVP);
    shader->setShaderParamFromMat4("M",M);
}

void GLWindow::setSelectedObject(RenderObject *object)
{
    m_selectedObject = object;
    m_selectedTransform.reset();
    if (m_selectedObject != NULL)
    {
        m_selectedTransform.m_transform = m_selectedObject->getTransform();
        m_selectedTransform.m_translation.set(m_selectedTransform.m_transform.m_30, m_selectedTransform.m_transform.m_31, m_selectedTransform.m_transform.m_32);
    }
}

void GLWindow::drawStrand(const Strand &strand)
{
    if (m_strandVAO != NULL)
    {
        m_strandVAO->bind();

        m_strandVAO->updateIndexedData(strand.m_ppos.size() * sizeof(ngl::Vec4),
                                       strand.m_ppos[0].m_x);

        m_strandVAO->setVertexAttributePointer(0, 4, GL_FLOAT, 0, 0);

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


    m_strandVAO->setIndexedData(strand.m_ppos.size() * sizeof(ngl::Vec4),
                               strand.m_ppos[0].m_x,
                               indices.size(),
                               &indices[0],
                               GL_UNSIGNED_INT);

    m_strandVAO->setVertexAttributePointer(0, 4, GL_FLOAT, 0, 0);
    m_strandVAO->setNumIndices(indices.size());

    m_strandVAO->draw();
    m_strandVAO->unbind();
}

void GLWindow::drawHairStrand(const ElasticRod& strand)
{
    if (m_strandVAO != NULL)
    {
        m_strandVAO->bind();

        m_strandVAO->updateIndexedData(strand.m_ppos.size() * sizeof(ngl::Vec4),
                                       strand.m_ppos[0].m_x);

        m_strandVAO->setVertexAttributePointer(0, 4, GL_FLOAT, 0, 0);

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

    m_strandVAO->setIndexedData(strand.m_ppos.size() * sizeof(ngl::Vec4),
                               strand.m_ppos[0].m_x,
                               indices.size(),
                               &indices[0],
                               GL_UNSIGNED_INT);

    m_strandVAO->setVertexAttributePointer(0, 4, GL_FLOAT, 0, 0);
    m_strandVAO->setNumIndices(indices.size());

    m_strandVAO->draw();
    m_strandVAO->unbind();
}

//----------------------------------------------------------------------------------------------------------------------
//This virtual function is called whenever the widget needs to be painted.
// this is our main drawing routine
void GLWindow::paintGL()
{
    if ((m_rotate || m_zoom) && m_selectedObject != NULL)
    {
        m_selectedObject->setTransform(m_selectedTransform.getTransform());
    }

    // clear the screen and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // grab an instance of the shader manager
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();
    // get an instance of the VBO primitives for drawing
    ngl::VAOPrimitives *prim=ngl::VAOPrimitives::instance();

    m_transform.reset();

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

    // draw spring lines
//    glLineWidth(0.05);
//    shader->use("Colour");
//    shader->setShaderParam4f("Colour",0.4,0.3,0.0,1.0);
//    // load transform stack
//    m_transform.reset();
//    loadMatricesToHairShader();

//    const std::vector<Hair*>& hList = m_scene->getHairs();
//    typedef std::vector<Hair*>::const_iterator HIter;
//    typedef std::vector<Strand>::const_iterator SIter;
//    for (HIter hit = hList.begin(); hit != hList.end(); hit++)
//    {
//        Hair* hair = *hit;
//        for (SIter it = hair->m_strands.begin(); it != hair->m_strands.end(); it++)
//        {

//            drawStrand((*it));
//        }
//    }

    // draw spring lines
    glLineWidth(0.05);
    shader->use("Colour");
    shader->setShaderParam4f("Colour", 0.8, 0.8, 0.0, 1.0);
    // load transform stack
    m_transform.reset();
    loadMatricesToHairShader();

    const std::vector<ElasticRod*>& strands = m_scene->getStrands();
    typedef std::vector<ElasticRod*>::const_iterator SIter;
    for (SIter it = strands.begin(); it != strands.end(); ++it)
    {
        drawHairStrand(**it);
    }

    shader->use("Phong");
//    shader->setShaderParam4f("Colour",0.4,0.4,0.4,1.0);
    m_transform.reset();
    loadMatricesToShader();
    prim->draw("grid");

}





//----------------------------------------------------------------------------------------------------------------------
void GLWindow::mouseMoveEvent ( QMouseEvent * _event )
{
    TransformTool& ttool = (m_selectedObject == NULL)? m_cameraTransform : m_selectedTransform;

    // note the method buttons() is the button state when event was called
    // this is different from button() which is used to check which button was
    // pressed when the mousePress/Release event is generated
    if(m_rotate && _event->buttons() == Qt::LeftButton)
    {
        ngl::Real dx = _event->x() - m_mouseX;
        ngl::Real dy = _event->y() - m_mouseY;
        ttool.m_phi = fmod(ttool.m_phi + dx, 360.0);
        ttool.m_theta = fmod(ttool.m_theta + dy, 360.0);

        m_mouseX = _event->x();
        m_mouseY = _event->y();
        // re-draw GL
        update();
    }

    if(m_zoom && _event->buttons() == Qt::RightButton)
    {
        ngl::Real dx = 2 * ngl::Real(_event->x() - m_mouseX) / width();
        ngl::Real dy = 2 * ngl::Real(_event->y() - m_mouseY) / height();
        ngl::Vec4 dirz = (m_cam->getLook() - m_cam->getEye());

        if (m_selectedObject != NULL)
        {
            dirz = dirz * m_cameraTransform.getTransform().inverse();
            dirz.normalize();
            ngl::Vec4 dirx = utils::EY.cross(dirz);
            dirx.normalize();
            ngl::Vec4 diry = dirz.cross(dirx);
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
        case Qt::MidButton:
            m_pan = true;
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
        case Qt::MidButton:
            m_pan = false;
            break;

        default:
            break;
    }
}

void GLWindow::timerEvent( QTimerEvent *_event)
{
    m_scene->update(1.0 / utils::UPS);
    updateGL();
}

void GLWindow::startSimTimer()
{

    m_timer = startTimer(m_timerValue);
}

void GLWindow::stopSimTimer()
{
    killTimer(m_timer);
}


// ============================ utility functions ==========================

void buildVAOs(const std::vector<Mesh*>& meshList, std::vector<ngl::VertexArrayObject*>& o_VAOList)
{
    assert(o_VAOList.size() == 0);
    o_VAOList.resize(meshList.size(), NULL);

    unsigned i = 0;
    for (std::vector<Mesh*>::const_iterator it = meshList.begin(); it != meshList.end(); ++it)
    {
        if ((*it)->getNFaces() == 0)
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

void feedVAO(const Mesh& mesh, ngl::VertexArrayObject& o_vao)
{
    o_vao.bind();
    // in this case we are going to set our data as the vertices above
    // now we set the attribute pointer to be 0 (as this matches vertIn in our shader)
    o_vao.setIndexedData(mesh.m_vertices.size() * sizeof(ngl::Vec4),
                        mesh.m_vertices[0].m_x,
                        mesh.m_vindices.size(),
                        &mesh.m_vindices[0],
                        GL_UNSIGNED_INT);

    o_vao.setVertexAttributePointer(0, 4, GL_FLOAT, 0, 0);

    // now we set the attribute pointer to be 2 (as this matches normalIn in our shader)
    o_vao.setIndexedData(mesh.m_normals.size() * sizeof(ngl::Vec4),
                        mesh.m_normals[0].m_x,
                        mesh.m_nindices.size(),
                        &mesh.m_nindices[0],
                        GL_UNSIGNED_INT);
    o_vao.setVertexAttributePointer(2, 4, GL_FLOAT, 0, 0);
    o_vao.setNumIndices(mesh.m_vindices.size());

    // now unbind
    o_vao.unbind();
}
