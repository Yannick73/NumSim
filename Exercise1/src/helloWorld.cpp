#include <iostream>
#include "Settings/Settings.hpp"

int main() {

	Settings test_settings;
    std::cout << "Hello world, tau: " << test_settings.tau << std::endl;
    test_settings.printSettings();
}
