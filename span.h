#pragma once

#include <vector>

struct ValueSpan
{
	ValueSpan(float min, float max, float step = 1) {
		auto l = min;
		while (l < max) {
			m_Latitudes.push_back(l);
			l += step;
		}

		m_Latitudes.push_back(max); // push max
	}

	// allows ranged-based for loops to work.
	auto begin() { return m_Latitudes.begin(); }
	auto end() { return m_Latitudes.end(); }

private:
	std::vector<float> m_Latitudes;
};