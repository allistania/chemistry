#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <unordered_map>
#include <regex>

using namespace std;

// Structure for representing a chemical component
struct Component {
    string name; // The name of the component (e.g. "H2O")
    map<string, int> elements; // Elements and their number (e.g. {"H": 2, "O": 1})
    int charge; // Charge (if it is an ion)
};

// Elemental charge database
unordered_map<string, int> elementCharges = {
    {"H", 1},  // hydrogen: +1
    {"O", -2},   // oxygen: -2
    {"Li", 1},  // lithium: +1
    {"Be", 2},   // beryllium: +2
    {"B", 3},   // boron: +3
    {"C", 2},   // carbon: +2
    {"F", -1},   // fluorine: -1
    {"Na", 1},  // sodium: +1
    {"Mg", 2},   // magnesium: +2
    {"Al", 3},   // aluminium: +3
    {"Si", -4},   // silicon: -4
    {"Si", 4},   // silicon: +4
    {"K", 1},   // potassium: +1
    {"Ca", 2},   // calcium: +2
    {"Sc", 3},   // scandium: +3
    {"Cl", -1},   // chlorine: -1
    {"Ti", -4},   // titanium: -4
    {"V", 5},   // vanadium: +5
    {"Cr", 6},   // chrome: +6
    {"Mn", 4},   // manganese: +4
    {"Fe", 2},   // iron: +2
    {"Co", 3},   // cobalt: +3
    {"Ni", 2},   // nickel: +2
    {"Cu", 2},   // copper: +2
    {"Zn", 2},   // zinc: +2
    {"Ce", 1},  // cesium: +1
    {"Ba", 2}   // barium: +2
    // need to add other elements
};

// Function to parse a component string and get the elements and their number
vector<pair<string, int>> parseElementCounts(const string& componentStr) {
    vector<pair<string, int>> elementCounts;
    regex elementRegex("([A-Z][a-z]*)(\\d*)"); // Regular expression for finding elements and their number

    sregex_iterator it(componentStr.begin(), componentStr.end(), elementRegex);
    sregex_iterator end;

    while (it != end) {
        string element = (*it)[1].str(); // element
        string countStr = (*it)[2].str(); // quantity

        int count = 1; // quantity
        if (!countStr.empty()) {
            count = stoi(countStr); // convert string to a number
        }

        elementCounts.push_back(make_pair(element, count));
        ++it;
    }

    return elementCounts;
}

// Function for parsing a chemical component
Component parseComponent(const string& componentStr) {
    Component component;
    component.name = componentStr;
    component.charge = 0;

    vector<pair<string, int>> elementCounts = parseElementCounts(componentStr);

    for (const auto& elementCount : elementCounts) {
        string element = elementCount.first;
        int count = elementCount.second;

        // Getting the charge of an element from a table
        if (elementCharges.count(element)) {
            component.elements[element] = count;
            component.charge += elementCharges[element] * count;
        } else {
            cerr << "Error: Unknown element " << element << endl;
        }
    }

    return component;
}

// Function to separate components of an equation
vector<Component> splitComponents(const string& side) {
    vector<Component> components;
    stringstream ss(side);
    string componentStr;
    while (getline(ss, componentStr, '+')) {
        Component component = parseComponent(componentStr);
        components.push_back(component);
    }
    return components;
}

// Function for outputting information about the component
void printComponent(const Component& component) {
    cout << component.name << " (";
    for (const auto& element : component.elements) {
        cout << element.first << ":" << element.second << " ";
    }
    cout << "Charge: " << component.charge << ")";
}

// Function for outputting information about the reactans
void printReactants(const Component& component) {
    cout << component.name << " (";
    for (const auto& element : component.elements) {
        cout << element.first << ":" << element.second << " ";
    }
    cout << "Charge: " << 0 << ")";
}


int main() {
    string equation;
    cout << "Enter a chemical equation (e.g., H2 + O2 = H2O): ";
    getline(cin, equation);

    // Separation of the equation into reactants and products
    size_t equalPos = equation.find('=');
    if (equalPos == string::npos) {
        cerr << "Invalid equation format." << endl;
        return 1;
    }

    string reactantsStr = equation.substr(0, equalPos);
    string productsStr = equation.substr(equalPos + 1);

    vector<Component> reactants = splitComponents(reactantsStr);
    vector<Component> products = splitComponents(productsStr);

// Output the equation
    cout << "Equation: " << equation << endl;

    // Displaying information about components
    cout << "Reactants: " << endl;
    for (const auto& component : reactants) {
        printReactants(component);
        cout << endl;
    }

    cout << "Products: " << endl;
    for (const auto& component : products) {
        printComponent(component);
        cout << endl;
    }

    return 0;
}
