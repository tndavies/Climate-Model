#include <datapack.h>
#include <iostream>
#include <span.h>
#include <cmath>
#include <sim.h>

int main(int argc, char* argv[])
{

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

	const float Year = 3.154e7; // 1 year in seconds
	const size_t MonthlySampleRate = 1e3;
	
	ValueSpan latitudes(-90, 90, 10);
	size_t pack_id = 1;
	
	for (const auto& lat : latitudes) 
	{
		Datapack pack;
		
		ValueSpan times(0, Year, Year / (12.0f * MonthlySampleRate));
		for (const auto& t : times) {

			const float flux = Sim::Calc_DiurnalSolarFlux(lat, Sim::Calc_Declination(t));
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