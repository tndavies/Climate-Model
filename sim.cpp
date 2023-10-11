#include <sim.h>
#include <cmath>

#define Radians(x) (0.01745329252f * x) // convert x from degrees to radians.

float Sim::Calc_Declination(float t)
{
	// Values for Earth ..
	const float Obliquity = Radians(23.45f); // sign controls which way we orbit sun.
	const float OrbitalAngVel = 1.98226e-7f;
	//

	const float Offset = 1.570796327f; // pi/2
	const float sine_decl = -1.0f * std::sinf(Obliquity) * std::cosf(OrbitalAngVel * t + Offset);

	return std::asinf(sine_decl);
}

float Sim::Calc_DiurnalSolarFlux(float latitude, float declination)
{
	const float pi = 3.14159265358979f;
	const float prefactor = 1360.0f / pi;
	const float lat_Rad = Radians(latitude);

	float H = 0.0f;
	const float tan_product = std::tanf(lat_Rad) * std::tanf(declination);
	if (tan_product > 1.0f) {
		H = pi;
	}
	else if (tan_product < -1.0f) {
		return 0.0f;
	}
	else {
		H = std::acosf(-tan_product);
	}
	
	float term0 = H * std::sinf(lat_Rad) * std::sinf(declination);
	float term1 = std::cosf(lat_Rad) * std::cosf(declination) * std::sinf(H);
	float diurnal_flux = prefactor * (term0 + term1);

	return(diurnal_flux);
}