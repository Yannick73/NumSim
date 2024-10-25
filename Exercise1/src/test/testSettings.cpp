#include <iostream>
#include "Settings/Settings.hpp"

int main() {

	Settings test_settings;
    std::cout << "Test settings loading procedure. Settings before loading file\n\n";
    test_settings.printSettings();
    test_settings.loadFromFile("../../lid_driven_cavity.txt");
    std::cout << "\nSettings after loading file\n\n";
    test_settings.printSettings();
}
