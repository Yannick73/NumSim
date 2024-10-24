#include "Settings.hpp"

#define PRINT_ARRAY(array) \
	std::cout << #array << ": \t" << array2str(array) << std::endl
#define PRINT_ARRAY_UNIT(array, unit) \
	std::cout << #array << ": \t" << array2str(array) << "\t[" << \
	unit << "]" << std::endl
#define PRINT_PARAM(param) \
	std::cout << #param << ": \t" << param <<  std::endl
#define PRINT_PARAM_UNIT(param, unit) \
	std::cout << #param << ": \t" << param << "\t" << unit <<  std::endl

#define LOAD_INT(lines, param) { int line_num = find_param(lines, param); \
	if(line_num > 0) { #param = extract_int(lines[line_num]); }}
#define LOAD_FLOAT(param) { int line_num = find_param(lines, param); \
	if(line_num > 0) { #param = atof(lines[line_num]); }}
#define LOAD_INT_ARRAY(param, dir, ind) \
	{ int line_num = find_param(lines, param##dir); \
	if(line_num > 0) { #param[ind] = extract_int(lines[line_num]); }}
#define LOAD_FLOAT_ARRAY(param, dir, ind) \
	{ int line_num = find_param(lines, param##dir); \
	if(line_num > 0) { #param[ind] = atof(lines[line_num]); }}
	
int extract_int(const std::string &input)
{
	// per convention starts each line with the param name followed by
	// "=" and then the value
	std::string value = input;
	
}

template <typename T, std::size_t N>
std::string array2str(const std::array<T, N> &input)
{
	std::stringstream stream;
	stream << "(";
	std::size_t last_elem = N - 1;
	
	for(int i = 0; i < N - 1; i++)
	{
		stream << input[i] << ", ";
	}
	stream << input[N - 1] << ")";
	
	return stream.str();
}

void Settings::loadFromFile(std::string filename)
{
	// load settings file
	std::ifstream settings_input(filename);
	// check, whether the file is valid and if not, throw an error
	if(!settings_input)
	{
		std::stringstream error_out;
		error_out << "Error: The file " << filename << \
			" could not be opened.";
		throw std::runtime_error(error_out.str());
	}
	
	// then load the file a line at a time
	std::vector<std::string> lines;
	std::string current_line;
	while(std::getline(settings_input, current_line))
	{
		lines.push_back(current_line);
	}
	
	LOAD_INT_ARRAY(nCells, X, 0);
	LOAD_INT_ARRAY(nCells, Y, 1);
	LOAD_FLOAT_ARRAY(physicalSize, X, 0);
	LOAD_FLOAT_ARRAY(physicalSize, Y, 1);
	LOAD_FLOAT(re);
	LOAD_FLOAT(endTime);
	LOAD_FLOAT(tau);
	LOAD_FLOAT(maximumDt);
	
	LOAD_FLOAT_ARRAY(g, X, 0);
	LOAD_FLOAT_ARRAY(g, Y, 1);
	
	// load bool
	LOAD_FLOAT(alpha);
	
	LOAD_FLOAT_ARRAY(dirichletBottom, X, 0);
	LOAD_FLOAT_ARRAY(dirichletBottom, Y, 1);
	LOAD_FLOAT_ARRAY(dirichletTop, X, 0);
	LOAD_FLOAT_ARRAY(dirichletTop, Y, 1);
	LOAD_FLOAT_ARRAY(dirichletLeft, X, 0);
	LOAD_FLOAT_ARRAY(dirichletLeft, Y, 1);
	LOAD_FLOAT_ARRAY(dirichletRight, X, 0);
	LOAD_FLOAT_ARRAY(dirichletRight, Y, 1);
	
	// load string
	LOAD_FLOAT(omega);
	LOAD_FLOAT(epsilon);
	LOAD_INT(maximumNumberOfIterations);
}

void Settings::printSettings()
{
	std::cout << "Settings used in simulation:" << std::endl;
	PRINT_ARRAY(nCells);
	PRINT_ARRAY_UNIT(physicalSize, "m");
	PRINT_PARAM(re);
	PRINT_PARAM_UNIT(endTime, "s");
	PRINT_PARAM(tau);
	PRINT_PARAM_UNIT(maximumDt, "s");
	
	PRINT_ARRAY_UNIT(g, "N");
	
	PRINT_PARAM(useDonorCell);
	PRINT_PARAM(alpha);
	
	PRINT_ARRAY_UNIT(dirichletBottom, "m/s");
	PRINT_ARRAY_UNIT(dirichletTop, "m/s");
	PRINT_ARRAY_UNIT(dirichletLeft, "m/s");
	PRINT_ARRAY_UNIT(dirichletRight, "m/s");
	
	PRINT_PARAM(pressureSolver);
	PRINT_PARAM(omega);
	PRINT_PARAM(epsilon);
	PRINT_PARAM(maximumNumberOfIterations);
}
