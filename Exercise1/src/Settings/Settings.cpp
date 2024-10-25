#include "Settings.hpp"

// print macro, which prints parameter name, value and unit (if applicable)
#define PRINT_ARRAY(array) \
	std::cout << #array << ": \t" << array2str(array) << std::endl
#define PRINT_ARRAY_UNIT(array, unit) \
	std::cout << #array << ": \t" << array2str(array) << "\t[" << \
	unit << "]" << std::endl
#define PRINT_PARAM(param) \
	std::cout << #param << ": \t" << param <<  std::endl
#define PRINT_PARAM_UNIT(param, unit) \
	std::cout << #param << ": \t" << param << "\t" << unit <<  std::endl

// wraps the param to a name and a string, to make user error less likely
#define LOAD_DOUBLE(param) load_double(&param, #param)
// same as load-double, but dir is the direction literal (aka 'X' or 'Y') used in the file,
// whereas ind is the index of the array/vectir
#define LOAD_DOUBLE_ARRAY(param, dir, ind) load_double(&(param[ind]), #param #dir)
#define LOAD_INTEGER(param) load_integer(&param, #param)
#define LOAD_INTEGER_ARRAY(param, dir, ind) load_integer(&(param[ind]), #param #dir)
#define LOAD_BOOLEAN(param) load_boolean(&param, #param)
#define LOAD_STRING(param) load_string(&param, #param)

std::string strclean(const std::string input)
{
	std::string output = input;
	// transform each char to lower case
	std::transform(output.begin(), output.end(), output.begin(),
    [](unsigned char c){ return std::tolower(c); });
	// remove all whitespace and tab characters
	output.erase(remove(output.begin(), output.end(), ' '), output.end());
	output.erase(remove(output.begin(), output.end(), 't'), output.end());

	return output;
}

std::string Settings::extract_string(const std::string name)
{
	std::string cleaned_name = strclean(name);
	std::string output = "";
	// iterate over each line
	for(std::string line_content : file_lines)
	{
		// if a line defines a parameter, it contains a '='
		std::size_t index = line_content.find_first_of('=');
		if(index != std::string::npos) 
		{
			std::string line_begin = line_content.substr(0, index);
			std::string line_end = line_content.substr(index+1);
			// if the basic string content is identical, return the basic value
			if(strclean(line_begin) == cleaned_name)
			{
				output = strclean(line_end);
				break;
			}
		}
	}

	return output;
}

void Settings::load_double(double *attribute, const std::string name)
{
	std::string value = extract_string(name);
	*attribute = atof(value.c_str());
}

void Settings::load_integer(int *attribute, const std::string name)
{
	std::string value = extract_string(name);
	*attribute = atoi(value.c_str());
}

void Settings::load_boolean(bool *attribute, const std::string name)
{
	std::string value = extract_string(name);
	std::string true_strings[] = {"true", "1", "yes"};
	bool value_true = false;
	// set the attribute to true, if it is equal to one of the possible true-values
	for(std::string cmp_str : true_strings)
	{
		if(value == cmp_str)
		{
			value_true = true;
			break;
		}
	}
	*attribute = value_true;
}

void Settings::load_string(std::string *attribute, const std::string name)
{
	*attribute = extract_string(name);
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
	
	std::string current_line;
	// iterate over all the lines
	while(std::getline(settings_input, current_line))
	{
		// remove any comments indicated by "#"
		std::size_t index = current_line.find_first_of('#');
		if(index != std::string::npos)
		{
			// only keep the string from beggining to the found '#'
			current_line = current_line.substr(0, index);
		}
		// remove all whitespaces and set to lowercase
		current_line = strclean(current_line);
		// only non-empty lines are relevant and are saved
		if(current_line != "")
		{
			file_lines.push_back(current_line);
		}
	}

	// load each relevant setting
	LOAD_INTEGER_ARRAY(nCells, X, 0);
	LOAD_INTEGER_ARRAY(nCells, Y, 1);
	LOAD_DOUBLE_ARRAY(physicalSize, X, 0);
	LOAD_DOUBLE_ARRAY(physicalSize, Y, 1);
	LOAD_DOUBLE(re);
	LOAD_DOUBLE(endTime);
	LOAD_DOUBLE(tau);
	LOAD_DOUBLE(maximumDt);
	
	LOAD_DOUBLE_ARRAY(g, X, 0);
	LOAD_DOUBLE_ARRAY(g, Y, 1);
	
	LOAD_BOOLEAN(useDonorCell);
	LOAD_DOUBLE(alpha);
	
	LOAD_DOUBLE_ARRAY(dirichletBottom, X, 0);
	LOAD_DOUBLE_ARRAY(dirichletBottom, Y, 1);
	LOAD_DOUBLE_ARRAY(dirichletTop, X, 0);
	LOAD_DOUBLE_ARRAY(dirichletTop, Y, 1);
	LOAD_DOUBLE_ARRAY(dirichletLeft, X, 0);
	LOAD_DOUBLE_ARRAY(dirichletLeft, Y, 1);
	LOAD_DOUBLE_ARRAY(dirichletRight, X, 0);
	LOAD_DOUBLE_ARRAY(dirichletRight, Y, 1);
	
	LOAD_STRING(pressureSolver);
	LOAD_DOUBLE(omega);
	LOAD_DOUBLE(epsilon);
	LOAD_INTEGER(maximumNumberOfIterations);
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
