#ifndef CONFIG_H
#define CONFIG_H

#include "ngl/Util.h"
#include "ngl/Vec4.h"

namespace config
{

    const ngl::Real sceneBoundingBoxSize = 60;
    const unsigned gridDivisions = 10;

    const unsigned nInitialObstacles = 5;
    const unsigned nInitialBoids = 1000;

    const ngl::Real boidInitialSpeedMin = 3;
    const ngl::Real boidMaxSpeed = 6;
    const ngl::Real boidMass = 2;
    const ngl::Real boidMaxTurningAngle = ngl::PI4;
    const ngl::Real boidPanicDist = 3;
    const ngl::Real boidNeighbourhoodDist = 6;
    const ngl::Real boidNeighbourhoodFOV = ngl::PI;
    const ngl::Real boidObstacleLookupMinDist = 8;

    const ngl::Real obstacleAvoidanceWeight = 1.0;
    const ngl::Real separationWeight = 0.7;
    const ngl::Real alignmentWeight = 0.8;
    const ngl::Real cohesionWeight = 0.7;
    const ngl::Real volumeConstraintWeight = 1.0;
    const ngl::Real wanderWeight = 1.0;

    extern bool g_useMultiThreading;
}

#endif // CONFIG_H
