#include <datapack.h>
#include <fstream>

void Datapack::add(float x, float y)
{
	m_Data.push_back({.x=x, .y=y});
}

void Datapack::dump(std::string datapack_name, std::string comment)
{
	const char CommentSymbol = '#';

	std::ofstream datapack;
	datapack.open(datapack_name.c_str(), std::ios::out | std::ios::trunc);

	if (comment.length()) {
		datapack << CommentSymbol << " ";
		datapack << comment;
		datapack << std::endl;
	}

	for (const auto& dp : m_Data) {
		datapack << dp.x << ',' << dp.y << std::endl;
	}

	datapack.close();
}