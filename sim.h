#pragma once

namespace Sim {
	float Calc_Declination(float t);
	float Calc_DiurnalSolarFlux(float latitude, float declination);
};