#ifndef __GL_WINDOW_H__
#define __GL_WINDOW_H__

#include "Scene.h"
#include "GLUtils.h"
#include "TransformTool.h"

#include <QOpenGLWidget>
#include <QString>


//----------------------------------------------------------------------------------------------------------------------
/// @class GLWindow
/// @brief our main gl window widget - all viewport drawing happens here
//----------------------------------------------------------------------------------------------------------------------
class GLViewport : public QOpenGLWidget
{

// must include this if we want to use Qt signals/slots
Q_OBJECT

public :
	explicit GLViewport(QWidget* parent = Q_NULLPTR, Qt::WindowFlags f = Qt::WindowFlags());
	~GLViewport();

/// @brief set the scene to render
/// @param[in] scene - the scene
	inline void setScene(Scene* scene)
	{
		m_scene = scene;
	}

	inline void setSelection(bool selection)
	{
		m_selection = selection;
//		m_selectionTransform.reset();
	}

	inline bool getSelection() const
	{
		return m_selection;
	}

	inline void setSelectionTransform(const mg::Matrix4D& transform)
	{
//		m_selectionTransform.reset();
//		m_selectionTransform.setTransform(transform);
//		updateGL();
	}

	inline mg::Matrix4D getSelectionTransform()
	{
		static mg::Matrix4D tm;
		return tm;
//		return m_selectionTransform.getTransform();
	}

protected:
/// @brief  The following methods must be implimented in the sub class
/// this is called when the window is created
	virtual void initializeGL();

/// @brief this is called whenever the window is re-sized
/// @param[in] w the width of the resized window
/// @param[in] h the height of the resized window
	virtual void resizeGL(int w, int h);

/// @brief this is the main gl drawing routine which is called whenever
/// the window needs to be re-drawn
	virtual void paintGL();

	virtual void mousePressEvent(QMouseEvent* event);
	virtual void mouseMoveEvent(QMouseEvent* event);
	virtual void mouseReleaseEvent(QMouseEvent* event);

private :
	typedef std::vector< std::unique_ptr< GLDrawable > > DrawList;

// ============= various draw related utility functions ===============
	void buildVAOs(const std::vector<Mesh*>& meshList, DrawList& o_drawList) const;
	void drawHairStrand(const ElasticRod& strand);

private:
/// @brief viewport camera
	Camera m_cam;
/// @brief reference grid
	GLDrawable m_refGrid;
/// @brief debug info overlay
//	TODO fix: ngl::Text* m_text;

/// @brief scene to render
	const Scene* m_scene = nullptr;
/// @brief list of drawable objects created for each mesh in scene
	DrawList m_drawList;

/// mouse controls
/// @brief flag that marks that transformations calculated on mouse movement apply to an object and not the camera
	bool m_selection = false;
/// @brief transform handle
	TransformHandle m_transformHdl;
/// @brief the previous x mouse value
	int m_mouseX = 0;
/// @brief the previous y mouse value
	int m_mouseY = 0;
};

#endif
