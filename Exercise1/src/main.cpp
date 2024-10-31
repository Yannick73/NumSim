#include <iostream>
#include <exception>
#include "Settings/Settings.hpp"

int main(int argc, char* argv[])
{
    Settings sim_settings;

    if(argc == 1) 
    {
        // argument 0 is the path to the binary
        std::cout << "No settings file specified, start sim with standard settings\n";
    }
    else if(argc == 2)
    {
        // argument 1 is the first optional argument
        std::string settings_path = argv[1];
        try {
            sim_settings.loadFromFile(settings_path);
            std::cout << "Loaded settings: " << settings_path << std::endl;
        } catch(std::runtime_error &e) {
            std::cerr << "Error loading file: " << e.what() << "\nStop program.\n";
            return 1;
        }
    }
    else
    {
        std::cerr << "The only possible argument is the path to a settings file. Stop program.\n";
        return 1;
    }

    sim_settings.printSettings();

    // do stuff

    return 0;
}