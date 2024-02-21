#ifndef GMGBOUNDARYTRAIT_H
#define GMGBOUNDARYTRAIT_H
#pragma once
#include "FirstOrderUpdate.h"
#include "FirstOrderKFVSBoundaryPackage.h"

template<typename FluxSolverType>
struct FluxSolverToBoundary
{

};

template<>
struct FluxSolverToBoundary<FirstOrderUpdate<NS, 1>>
{
	using BoundaryType = typename FirstOrderKFVSBoundaryPackage;
};

template<>
struct FluxSolverToBoundary<FirstOrderUpdate<NS, 2>>
{
	using BoundaryType = typename FirstOrderKFVSBoundaryPackage;
};

#endif;