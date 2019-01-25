#pragma once

#include "GLWindow.h"
#include "AnimationBuffer.h"
#include <QMainWindow>
#include <QTimer>

namespace Ui
{
class MainWindow;
}

class MainWindow : public QMainWindow
{
	Q_OBJECT
public:
	explicit MainWindow( QWidget *parent = nullptr );
	~MainWindow();

	void setScene( Scene& scene );

private slots :
	/// @brief slot to toggle the sim on and off
	void toggleSim(bool s);

	/// @brief timer event trigered by updateTimer and forward btn
	void updateEvent();

	/// @brief timer event trigered by recordTimer
	void recordEvent();

	/// @brief slot to toggle the recording on and off
	void toggleRecord(bool s);

	/// @brief set update timer interval
	void setTimerUpdateDuration(int ms);

	/// @brief slot to pick the selected object
	void selectRenderObject(int index);

	/// @brief slot to export current state of simulation
	void exportSim();

	void selectExportDirectory();

	void setBendingStiffness(double val);
	void setTwistingStiffness(double val);
	void setMaxElasticForce(double val);
	void setDrag(double val);
	void setPBDIter(int val);

	void selectMinimizationMethod(int idx);
	void setMinimizationTolerance(double val);
	void setMinimizationMaxIter(int val);

	void toggleCollisions(bool val);
	void toggleSelfInterations(bool val);
	void setSelfStiction(double val);
	void setSelfRepusion(double val);

private:
	void updateUI();

	void exportFrame(unsigned frame);

private:
	std::unique_ptr<Ui::MainWindow> m_ui; ///< UI layout
	GLViewport* m_gl = nullptr; ///< GL widget
	Scene* m_scene = nullptr; ///< scene to render
	QTimer m_updateTimer; ///< update timer
	/// @brief selected object to which transformations calculated on
	/// mouse movement are applied if nullptr camera is modified
	RenderObject* m_selectedObject = nullptr;

	bool m_exporting;
	unsigned m_frame;

	QString m_exportDir;
	QString m_exportGeoDir;
	QString m_exportCurvDir;

	QTimer m_recordTimer;
	AnimationBuffer m_animationBuffer;
};
