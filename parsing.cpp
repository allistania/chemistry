#include "functions.h"

vector<SavedReaction> savedReactions;
// Elements with fixed oxidation states
unordered_map<string, int> fixedOxidationStates = {
    {"H", 1}, {"O", -2}, {"F", -1}, {"Li", 1}, {"Be", 2}, {"B", 3},
    {"Na", 1}, {"Mg", 2}, {"Al", 3}, {"K", 1}, {"Ca", 2}, {"Sc", 3},
    {"Cl", -1}, {"Ti", 4}, {"V", 5}, {"Cr", 6}, {"Zn", 2}, {"Ce", 1},
    {"Ba", 2}, {"Br", -1}, {"I", -1}, {"Ag", 1}, {"Au", 3}
};

// Elements with variable oxidation states
unordered_map<string, vector<int>> variableOxidationStates = {
    {"C", {2, 4}}, {"Si", {-4, 4}}, {"N", {-3, 3, 5}}, {"P", {-3, 3, 5}},
    {"S", {-2, 4, 6}}, {"Mn", {2, 4, 6, 7}}, {"Fe", {2, 3}}, {"Co", {2, 3}},
    {"Ni", {2}}, {"Cu", {1, 2}}, {"Pb", {2, 4}}, {"Sn", {2, 4}}
};

// Function to calculate the oxidation state of an element in a compound
int calculateOxidationState(const string& element, const map<string, int>& elements, int totalCharge) {
    int knownChargeSum = 0;
    for (const auto& el : elements) {
        if (el.first != element) {
            if (fixedOxidationStates.count(el.first)) {
                knownChargeSum += fixedOxidationStates[el.first] * el.second;
            } else if (variableOxidationStates.count(el.first)) {
                knownChargeSum += variableOxidationStates[el.first][0] * el.second;
            }
        }
    }
    // Calculate the oxidation state for the current element
    return (totalCharge - knownChargeSum) / elements.at(element);
}

// Function to calculate the organiic oxidation state of an element in a compound
int calculateOrganicOxidationState(const string& element, const map<string, int>& elements, int totalCharge) {
    if (element == "C") {
        int knownChargeSum = 0;
        for (const auto& el : elements) {
            if (el.first != "C") {
                if (el.first == "H") knownChargeSum += 1 * el.second;
                else if (el.first == "O") knownChargeSum += -2 * el.second;
                else if (fixedOxidationStates.count(el.first)) {
                    knownChargeSum += fixedOxidationStates[el.first] * el.second;
                }
            }
        }
        return (totalCharge - knownChargeSum) / elements.at("C");
    }
    return 0;
}

vector<pair<string, int>> parseElementCounts(const string& componentStr) {
    vector<pair<string, int>> elementCounts;
    regex elementRegex("([A-Z][a-z]*)(\\d*)|(\\(([^)]+)\\))(\\d*)");

    sregex_iterator it(componentStr.begin(), componentStr.end(), elementRegex);
    sregex_iterator end;

    while (it != end) {
        smatch match = *it;
        if (match[1].matched) { // Regular element (no parentheses)
            string element = match[1].str();
            string countStr = match[2].str();
            int count = countStr.empty() ? 1 : stoi(countStr);
            elementCounts.push_back({element, count});
        }
        else if (match[3].matched) { // Group in parentheses
            string groupContent = match[4].str();
            string groupCountStr = match[5].str();
            int groupCount = groupCountStr.empty() ? 1 : stoi(groupCountStr);
            
            // Recursively parse the content within the parentheses
            auto innerElements = parseElementCounts(groupContent);
            for (auto& el : innerElements) {
                el.second *= groupCount; // Multiply by the group coefficient
                elementCounts.push_back(el);
            }
        }
        ++it;
    }
    return elementCounts;
}

// Function for parsing a chemical component
Component parseComponent(const string& componentStr) {
    Component component;
    component.name = componentStr;
    component.charge = 0; // Assume neutral charge unless specified otherwise
    component.coefficient = 1.0; // Default coefficient

    auto elementCounts = parseElementCounts(componentStr);
    for (const auto& ec : elementCounts) {
        component.elements[ec.first] += ec.second; // Use += to sum identical elements
    }

    if (component.elements.count("C") && component.elements["C"] > 0) {
        component.isOrganic = true;
        
        for (const auto& el : component.elements) {
            const string& elementName = el.first;
            if (elementName == "C") {
                component.oxidationStates["C"] = calculateOrganicOxidationState("C", component.elements, component.charge);
            }
            else if (elementName == "H") {
                component.oxidationStates["H"] = 1;
            }
            else if (elementName == "O") {
                component.oxidationStates["O"] = -2;
            }
            else {
                // Handle other elements (fixed or variable)
                if (fixedOxidationStates.count(elementName)) {
                    component.oxidationStates[elementName] = fixedOxidationStates[elementName];
                } else if (variableOxidationStates.count(elementName)) {
                    component.oxidationStates[elementName] = calculateOxidationState(elementName, component.elements, component.charge);
                }
            }
        }
    }
    else {
        for (const auto& el : component.elements) {
            const string& elementName = el.first;
            if (fixedOxidationStates.count(elementName)) {
                component.oxidationStates[elementName] = fixedOxidationStates[elementName];
            } else if (variableOxidationStates.count(elementName)) {
                component.oxidationStates[elementName] = calculateOxidationState(elementName, component.elements, component.charge);
            }
        }
    }
    return component;
}

vector<Component> splitComponents(const string& side) {
    vector<Component> components;
    stringstream ss(side);
    string componentStr;
    while (getline(ss, componentStr, '+')) {
        componentStr.erase(0, componentStr.find_first_not_of(" \t"));
        componentStr.erase(componentStr.find_last_not_of(" \t") + 1);
        if (!componentStr.empty()) {
            components.push_back(parseComponent(componentStr));
        }
    }
    return components;
}

// Function to collect all unique elements from reactants and products
vector<string> collectAllElements(const vector<Component>& reactants, const vector<Component>& products) {
    map<string, bool> elementsMap;
    for (const auto& c : reactants) {
        for (const auto& el : c.elements) elementsMap[el.first] = true;
    }
    for (const auto& c : products) {
        for (const auto& el : c.elements) elementsMap[el.first] = true;
    }
    vector<string> elements;
    for (const auto& pair : elementsMap) elements.push_back(pair.first);
    return elements;
}