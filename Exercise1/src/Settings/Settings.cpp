#include "Settings.hpp"

template <typename T, std::size_t N>
std::string array2str(const std::array<T, N> &input)
{
	std::stringstream stream;
	stream << "[";
	std::size_t last_elem = N - 1;
	
	for(int i = 0; i < N - 1; i++)
	{
		stream << input[i] << ", ";
	}
	stream << input[N - 1] << "]";
	
	return stream.str();
}

void Settings::loadFromFile(std::string filename)
{
	
}

void Settings::printSettings()
{
	std::cout << "nCells: " << array2str(nCells) << std::endl;
}
