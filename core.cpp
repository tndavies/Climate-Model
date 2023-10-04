#include <datapack.h>
#include <iostream>
#include <span.h>
#include <cmath>

float CalcDeclinationAngle(float obliquity, float angular_velocity, float t)
{
	const float PhaseShift = 1.57079632679490f; // pi/2.
	
	float sine_decl = -std::sinf(obliquity) * std::cosf(angular_velocity * t + PhaseShift);

	return std::asinf(sine_decl);
}

float CalcHalfDayLength(float latitude_band, float declination)
{
	float cosine_H = -std::tanf(latitude_band) * std::tanf(declination);

	return std::acosf(cosine_H);
}

float CalcAveragedSolarFlux(float declination, float latitude_band, float half_day_length)
{
	const float prefactor = 432.90144521f; // (flux @ 1 Au) / pi.
	
	const float term0 = half_day_length * std::sinf(latitude_band) * std::sinf(declination);
	const float term1 = std::cosf(latitude_band) * std::cosf(declination) * std::sinf(half_day_length);

	return prefactor * (term0 + term1);
}

int main(int argc, char* argv[])
{
	// Earth's orbital parameters.
	const float Obliquity = 23.5f;
	const float Omega = 1.983e-10f; 

	ValueSpan lats(-90, 90, 10);
	ValueSpan times(0, 3e8, 1e5);

	Datapack flux_data;

	for (const auto& lat : lats) 
	{
		for (const auto& t : times) {

			float decl = CalcDeclinationAngle(Obliquity, Omega, t);
			float half_day = CalcHalfDayLength(lat, decl);
			float S = CalcAveragedSolarFlux(decl, lat, half_day);

			flux_data.add(t, S);
		}
	}

	flux_data.dump("test.dp", "time (x) vs Diurnally averaged incident solar flux (y).");

	return(0);
}