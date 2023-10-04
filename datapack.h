#pragma once

#include <vector>
#include <string>

struct Datapack
{
	void add(float x, float y);
	void dump(std::string datapack_name, std::string comment = "");

private:
	struct DataPoint {
		float x, y;
	};

	std::vector<DataPoint> m_Data;
};