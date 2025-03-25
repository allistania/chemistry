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
    map<string, int> oxidationStates; // Oxidation states of elements
    int charge; // Charge (if it is an ion)
};

// Elements with fixed oxidation states
unordered_map<string, int> fixedOxidationStates = {
    {"H", 1},   // hydrogen: +1
    {"O", -2},  // oxygen: -2
    {"F", -1},  // fluorine: -1
    {"Li", 1},  // lithium: +1
    {"Be", 2},  // beryllium: +2
    {"B", 3},   // boron: +3
    {"Na", 1},  // sodium: +1
    {"Mg", 2},  // magnesium: +2
    {"Al", 3},  // aluminium: +3
    {"K", 1},   // potassium: +1
    {"Ca", 2},  // calcium: +2
    {"Sc", 3},  // scandium: +3
    {"Cl", -1}, // chlorine: -1
    {"Ti", 4},  // titanium: +4
    {"V", 5},   // vanadium: +5
    {"Cr", 6},  // chrome: +6
    {"Zn", 2},  // zinc: +2
    {"Ce", 1},  // cesium: +1
    {"Ba", 2},  // barium: +2
    {"Br", -1}, // bromine: -1
    {"I", -1},  // iodine: -1
    {"Ag", 1},  // silver: +1
    {"Au", 3}   // gold: +3
};

// Elements with variable oxidation states
unordered_map<string, vector<int>> variableOxidationStates = {
    {"C", {2, 4}},   // carbon: +2, +4
    {"Si", {-4, 4}}, // silicon: -4, +4
    {"N", {-3, 3, 5}}, // nitrogen: -3, +3, +5
    {"P", {-3, 3, 5}}, // phosphorus: -3, +3, +5
    {"S", {-2, 4, 6}}, // sulfur: -2, +4, +6
    {"Mn", {2, 4, 6, 7}}, // manganese: +2, +4, +6, +7
    {"Fe", {2, 3}},   // iron: +2, +3
    {"Co", {2, 3}},   // cobalt: +2, +3
    {"Ni", {2}},      // nickel: +2
    {"Cu", {1, 2}},   // copper: +1, +2
    {"Pb", {2, 4}},   // lead: +2, +4
    {"Sn", {2, 4}}    // tin: +2, +4
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

// Function to calculate the oxidation state of an element in a compound
int calculateOxidationState(const string& element, const map<string, int>& elements, int totalCharge) {
    int knownChargeSum = 0;
    for (const auto& el : elements) {
        if (el.first != element) {
            if (fixedOxidationStates.count(el.first)) {
                knownChargeSum += fixedOxidationStates[el.first] * el.second;
            } else {
                // For variable oxidation states, assume the first value in the list
                knownChargeSum += variableOxidationStates[el.first][0] * el.second;
            }
        }
    }
    // Calculate the oxidation state for the current element
    return (totalCharge - knownChargeSum) / elements.at(element);
}

// Function for parsing a chemical component
Component parseComponent(const string& componentStr) {
    Component component;
    component.name = componentStr;
    component.charge = 0; // Assume neutral charge unless specified otherwise

    vector<pair<string, int>> elementCounts = parseElementCounts(componentStr);

    for (const auto& elementCount : elementCounts) {
        string element = elementCount.first;
        int count = elementCount.second;

        component.elements[element] = count;
    }

    // Calculate oxidation states for all elements
    for (const auto& el : component.elements) {
        if (fixedOxidationStates.count(el.first)) {
            // Use fixed oxidation state
            component.oxidationStates[el.first] = fixedOxidationStates[el.first];
        } else if (variableOxidationStates.count(el.first)) {
            // Calculate oxidation state for elements with variable states
            int oxidationState = calculateOxidationState(el.first, component.elements, component.charge);
            component.oxidationStates[el.first] = oxidationState;
        } else {
            cerr << "Error: Unknown element " << el.first << endl;
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
    cout << "  Component: " << component.name << endl;
    for (const auto& element : component.elements) {
        string elementName = element.first;
        int count = element.second;
        int oxidationState = component.oxidationStates.at(elementName);
        cout << "    Element: " << elementName << " (Charge: " << oxidationState << ", Atoms: " << count << ")" << endl;
    }
    cout << "    Total charge of the component: " << component.charge << endl;
}

int main() {
    string equation;
    cout << "Enter a chemical equation (e.g., H2 + O2 = H2O): ";
    getline(cin, equation);

    // Separation of the equation into reactants and products
    size_t equalPos = equation.find('=');
    if (equalPos == string::npos) {
        cerr << "Invalis equation format." << endl;
        return 1;
    }

    string reactantsStr = equation.substr(0, equalPos);
    string productsStr = equation.substr(equalPos + 1);

    vector<Component> reactants = splitComponents(reactantsStr);
    vector<Component> products = splitComponents(productsStr);

    // Output the equation
    cout << "Equation: " << equation << endl;

    // Displaying information about reactants
    cout << "Reactants: " << endl;
    for (const auto& component : reactants) {
        printComponent(component);
    }

    // Displaying information about products
    cout << "Products: " << endl;
    for (const auto& component : products) {
        printComponent(component);
    }

    return 0;
}
