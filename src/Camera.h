#pragma once

#include "config.h"

class Camera
{
public:
	enum class ProjectionType
	{
		Perspective = 0,
		Orthographic,
		Sentinel,
	};

	Camera() { resetDefaults(); }

	const mg::Vec3D& getEye() const { return m_eye; }
	const mg::Vec3D& getAt() const { return m_at; }
	const mg::Vec3D& getUp() const { return m_up; }
	const mg::Matrix4D& getPMatrix() const { return m_projection; }
	const mg::Matrix4D& getVMatrix() const { return m_view; }
	mg::Matrix4D getVPMatrix() const { return m_projection * m_view; }
	ProjectionType getProjectionType() const { return m_projectionType; }

	void lookAt(const mg::Vec3D& eye);
	void lookAt(const mg::Vec3D& eye, const mg::Vec3D& at);
	void lookAt(const mg::Vec3D& eye, const mg::Vec3D& at, const mg::Vec3D& up);
	void ortho(mg::Real left, mg::Real right, mg::Real bottom, mg::Real top, mg::Real near, mg::Real far);
	void perspective(mg::Real yfov, mg::Real aspect, mg::Real near, mg::Real far);

	void rotate(mg::Real dx, mg::Real dy);
	void pan(mg::Real dx, mg::Real dy);
	void zoom(mg::Real dx, mg::Real dy);
	void resetDefaults();

private:
	ProjectionType m_projectionType;
	mg::Matrix4D m_projection;
	mg::Matrix4D m_view;

	mg::Vec3D m_eye;
	mg::Vec3D m_at;
	mg::Vec3D m_up;
};
