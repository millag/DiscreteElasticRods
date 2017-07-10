#include "Camera.h"


void Camera::resetDefaults()
{
	lookAt( mg::Vec3D(0.f, 0.f, 1.f), mg::Vec3D(0.f, 0.f, 0.f), mg::Vec3D(0.f, 1.f, 0.f ));
	ortho(-1.f, 1.f, -1.f, 1.f, -1.f, 1.f);
}

void Camera::lookAt(const mg::Vec3D& eye)
{
	m_eye = eye;
	mg::matrix_look_at_RH(m_view, m_eye, m_at, m_up);
}

void Camera::lookAt(const mg::Vec3D& eye, const mg::Vec3D& at)
{
	m_eye = eye;
	m_at = at;
	mg::matrix_look_at_RH(m_view, m_eye, m_at, m_up);
}

void Camera::lookAt(const mg::Vec3D& eye, const mg::Vec3D& at, const mg::Vec3D& up)
{
	m_eye = eye;
	m_at = at;
	m_up = up;
	mg::matrix_look_at_RH(m_view, m_eye, m_at, m_up);
}

void Camera::ortho(mg::Real left, mg::Real right, mg::Real bottom, mg::Real top, mg::Real near, mg::Real far)
{
	m_projectionType = ProjectionType::Orthographic;
	mg::matrix_orthographic_RH(m_projection, left, right, bottom, top, near, far, mg::z_clip_zero);
}

void Camera::perspective(mg::Real yfov, mg::Real aspect, mg::Real near, mg::Real far)
{
	m_projectionType = ProjectionType::Perspective;
	mg::matrix_perspective_yfov_RH(m_projection, yfov, aspect, near, far, mg::z_clip_zero);
}

void  Camera::rotate(mg::Real dx, mg::Real dy)
{
	mg::Matrix4D rot;
	rot.identity();
	mg::matrix_rotate_about_local_y(rot, -dx * mg::Constants::two_pi());
	mg::matrix_rotate_about_local_x(rot, -dy * mg::Constants::two_pi());
	lookAt( m_at + mg::transform_vector(rot, m_eye - m_at));
}

void  Camera::pan(mg::Real dx, mg::Real dy)
{
	const auto dirz = m_at - m_eye;
	auto dirx = mg::cross(m_up, dirz);
	auto diry = mg::cross(dirz, dirx);
	dirx.normalize();
	diry.normalize();
	const auto offset = dx * dirx + dy * diry;
	lookAt( m_eye + offset, m_at + offset );
}

void  Camera::zoom(mg::Real dx, mg::Real dy)
{
	const auto dirz = m_eye - m_at;
	lookAt( m_eye + dy * dirz );
}
