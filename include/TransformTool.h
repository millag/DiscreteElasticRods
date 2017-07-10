#ifndef TRANSFORMTOOL_H
#define TRANSFORMTOOL_H

#include "Camera.h"


class TransformHandle
{
public:
	enum TransformMode
	{
		TM_None = 0,
		TM_Rotation,
		TM_Scale,
		TM_Translation,
		TM_Count
	};

	TransformHandle(Camera* cam = nullptr):
		m_cam(cam)
	{ }
	~TransformHandle()
	{ }

	bool isActive() const { return m_cam && m_mode != TM_None; }
	TransformMode getMode() const { return m_mode; }
	void setMode(TransformMode mode) { m_mode = mode; }

	void update(mg::Real dx, mg::Real dy)
	{
		if (!m_cam)
		{
			return;
		}
		switch (m_mode)
		{
			case TM_Rotation:
			{
				m_cam->rotate(dx, dy);
				break;
			}
			case TM_Translation:
			{
				m_cam->pan(dx, dy);
				break;
			}
			case TM_Scale:
			{
				m_cam->zoom(dx, dy);
				break;
			}
		}
	}

private:
/// @brief transformation mode
	TransformMode m_mode = TM_None;
/// @brief transform object
	Camera* m_cam = nullptr;
};

#endif // TRANSFORMTOOL_H
