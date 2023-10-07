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

	return std::asinf(sine_decl); // @speed: can we use small angle approx?
}

float calc_half_day(float latitude, float declination)
{
	const auto f1 = -std::tanf(Radians(latitude));
	const auto f2 = std::tanf(declination);
	const float cosine_hd = f1 * f2;

	return std::acosf(cosine_hd);
}

float calc_diurnal_flux(float latitude, float declination, float half_day)
{
	const float prefactor = 432.90144521f; // Q_Naught / pi.
	const float term0 = half_day * std::sinf(Radians(latitude)) * std::sinf(declination);
	const float term1 = std::cosf(Radians(latitude)) * std::cosf(declination) * std::sinf(half_day);

	return prefactor * (term0 + term1);
}

int main(int argc, char* argv[])
{
	const float Earth_Obliquity = 23.45f;
	const float Earth_AngVel = 1.98226e-7f;
	const float Year = 3e7; // 1 year in seconds

#ifdef DATAPACK_DECLINATION
	const size_t MonthlySampleRate = 3;
	ValueSpan times(0, Year, Year / (12.0f * MonthlySampleRate));

	Datapack pack;
	for (const auto& t : times) {
		float decl = calc_declination(Earth_Obliquity, Earth_AngVel, t);
		pack.add(t, decl);
	}
	pack.dump("datapacks/declination.dp", "declination angle vs. time");
#endif

	ValueSpan latitudes(-90, 90, 10);
	const size_t MonthlySampleRate = 1e2;
	size_t pack_id = 1;
	
	for (const auto& lat : latitudes) 
	{
		Datapack pack;
		
		ValueSpan times(0, Year, Year / (12.0f * MonthlySampleRate));
		for (const auto& t : times) {

			const float decl = calc_declination(Earth_Obliquity, Earth_AngVel, t);
			const float half_day = calc_half_day(lat, decl); // @fix: we have a Nan issue here!
			const float flux = calc_diurnal_flux(lat, decl, half_day);

			pack.add(t, flux);
		}
		
		std::string pack_name("datapacks/");
		pack_name += "latitude_pack_" + std::to_string(pack_id++) + ".dp";
		
		std::string pack_title("latitude: ");
		pack_title += std::to_string(lat);

		pack.dump(pack_name, pack_title);
	}

	return(0);
}