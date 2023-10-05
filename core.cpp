#include <datapack.h>
#include <iostream>
#include <span.h>
#include <cmath>

#define Radians(x) (0.01745329252f * x)

float calc_declination(float obliquity, float omega, float t)
{
	const float obliquity_term = -std::sinf( Radians(obliquity) );
	const float orbit_term = std::sinf(omega * t);
	const float sine_decl = obliquity_term * orbit_term;

	return sine_decl; // small angle approx.
}

int main(int argc, char* argv[])
{
	const float Earth_Obliquity = 23.45f;
	const float Earth_AngVel = 1.98226e-7f;
	const float Year = 3e7; // 1 year in seconds

	const size_t MonthlySampleRate = 3;
	ValueSpan times(0, Year, Year / (12.0f * MonthlySampleRate));

	Datapack pack;

	for (const auto& t : times) {
		float decl = calc_declination(Earth_Obliquity, Earth_AngVel, t);
		pack.add(t, decl);
	}

	pack.dump("datapacks/declination.dp", "declination angle vs. time");

	return(0);
}