#include <iostream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

// Function for parsing chemical equation
vector<string> parseEquation(const string& equation) {
    vector<string> components;
    stringstream ss(equation);
    string component;
    while (getline(ss, component, '=')) {
        components.push_back(component);
    }
    return components;
}

// Function to separate components of an equation
vector<string> splitComponents(const string& side) {
    vector<string> components;
    stringstream ss(side);
    string component;
    while (getline(ss, component, '+')) {
        components.push_back(component);
    }
    return components;
}

// Function for outputting components separated by commas
void printComponents(const vector<string>& components) {
    for (size_t i = 0; i < components.size(); ++i) {
        cout << components[i];
        if (i < components.size() - 1) {
            cout << ", ";
        }
    }
    cout << endl;
}

int main() {
    string equation;
    cout << "Enter a chemical equation (e.g., H2 + O2 = H2O): ";
    getline(cin, equation);

// Parse equation
    vector<string> sides = parseEquation(equation);
    if (sides.size() != 2) {
        cerr << "Invalid equation format." << endl;
        return 1;
    }

// Separation into reactants and products
    vector<string> reactants = splitComponents(sides[0]);
    vector<string> products = splitComponents(sides[1]);

// Output the equation
    cout << "Equation: " << equation << endl;
// Output the reactants and products
    cout << "Reactants: ";
    printComponents(reactants);
    cout << "Products: ";
    printComponents(products);

    return 0;
}