#include "Scene.h"

void Scene::initialize()
{
	for ( auto i = 0u; i < m_meshes.size(); ++i )
	{
		m_meshes[i]->setId( i );
	}

	for ( auto i = 0u; i < m_renderObjects.size(); ++i )
	{
		m_renderObjects[i]->setId( i );
	}

	for ( auto i = 0u; i < m_hairs.size(); ++i )
	{
		m_hairs[i]->setId( i );
	}
}


void Scene::update(mg::Real dt)
{
	for ( auto it = m_hairs.begin(); it != m_hairs.end(); ++it )
	{
		(*it)->update( dt );
	}

//    m_spiral->update(dt);
}
