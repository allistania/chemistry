#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <unordered_map>
#include <regex>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>

using namespace std;

// Structure for representing a chemical component
struct Component {
    string name;
    map<string, int> elements;
    map<string, int> oxidationStates;
    int charge;
    double coefficient;
    bool isOrganic;
};

// Structure for storing saved reactions
struct SavedReaction {
    string originalEquation;
    string balancedEquation;
    vector<Component> reactants;
    vector<Component> products;
};

// Global vector to store saved reactions
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

// Function to display help information
void displayHelp() {
    cout << "Chemical Equation Balancer - Help\n";
    cout << "Usage:\n";
    cout << "  ./program [options] [equation]\n\n";
    cout << "Options:\n";
    cout << "  --help       Show this help message\n";
    cout << "  -v           Verbose output (show detailed information)\n\n";
    cout << "Examples:\n";
    cout << "  ./program \"H2 + O2 = H2O\"\n";
    cout << "  ./program -v \"Fe + O2 = Fe2O3\"\n";
    cout << "  ./program (launches interactive mode)\n\n";
    cout << "Interactive mode commands:\n";
    cout << "  1 - Balance a new equation\n";
    cout << "  2 - View saved reactions\n";
    cout << "  3 - Delete a saved reaction\n";
    cout << "  4 - Exit program\n";
}

// Function to parse a component string and get the elements and their number
vector<pair<string, int>> parseElementCounts(const string& componentStr) {
    vector<pair<string, int>> elementCounts;
    regex elementRegex("([A-Z][a-z]*)(\\d*)");

    sregex_iterator it(componentStr.begin(), componentStr.end(), elementRegex);
    sregex_iterator end;

    while (it != end) {
        string element = (*it)[1].str();
        string countStr = (*it)[2].str();
        int count = countStr.empty() ? 1 : stoi(countStr);
        elementCounts.push_back({element, count});
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

// Function for parsing a chemical component
Component parseComponent(const string& componentStr) {
    Component component;
    component.name = componentStr;
    component.charge = 0; // Assume neutral charge unless specified otherwise
    component.coefficient = 1.0; // Default coefficient

    auto elementCounts = parseElementCounts(componentStr);
    for (const auto& ec : elementCounts) {
        component.elements[ec.first] = ec.second;
    }

    if (component.elements.count("C") && component.elements["C"] > 0) {
        component.isOrganic = true;
        
        for (const auto& el : component.elements) {
            if (el.first == "C") {
                component.oxidationStates["C"] = calculateOrganicOxidationState("C", component.elements, component.charge);
            }
            else if (el.first == "H") {
                component.oxidationStates["H"] = 1;
            }
            else if (el.first == "O") {
                component.oxidationStates["O"] = -2;
            }
            else if (fixedOxidationStates.count(el.first)) {
                component.oxidationStates[el.first] = fixedOxidationStates[el.first];
            }
        }
    }
    else {
        for (const auto& el : component.elements) {
            if (fixedOxidationStates.count(el.first)) {
                component.oxidationStates[el.first] = fixedOxidationStates[el.first];
            }
            else if (variableOxidationStates.count(el.first)) {
                component.oxidationStates[el.first] = calculateOxidationState(el.first, component.elements, component.charge);
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

// Function to build the matrix for the system of equations
vector<vector<double>> buildEquationMatrix(const vector<string>& elements,
                                         const vector<Component>& reactants,
                                         const vector<Component>& products) {
    bool isOrganicReaction = false;
    // Number of variables (coefficients) is reactants.size() + products.size()
    int numVars = reactants.size() + products.size();
    vector<vector<double>> matrix(elements.size(), vector<double>(numVars + 1, 0.0));

    // For reactants (left side of equation)
    for (int i = 0; i < reactants.size(); ++i) {
        for (int j = 0; j < elements.size(); ++j) {
            if (reactants[i].elements.count(elements[j])) {
                matrix[j][i] = reactants[i].elements.at(elements[j]);
            }
        }
    }

    // For products (right side of equation, with negative sign)
    for (int i = 0; i < products.size(); ++i) {
        for (int j = 0; j < elements.size(); ++j) {
            if (products[i].elements.count(elements[j])) {
                matrix[j][reactants.size() + i] = -products[i].elements.at(elements[j]);
            }
        }
    }
    if (isOrganicReaction) {
        vector<double> carbonOxidationRow(numVars + 1, 0.0);
        for (int i = 0; i < reactants.size(); ++i) {
            if (reactants[i].oxidationStates.count("C")) {
                carbonOxidationRow[i] = reactants[i].oxidationStates.at("C") * reactants[i].elements.at("C");
            }
        }
        for (int i = 0; i < products.size(); ++i) {
            if (products[i].oxidationStates.count("C")) {
                carbonOxidationRow[reactants.size() + i] = -products[i].oxidationStates.at("C") * products[i].elements.at("C");
            }
        }
        matrix.push_back(carbonOxidationRow);
    }

    return matrix;
}

// Function to perform Gaussian elimination
void gaussianElimination(vector<vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    for (int i = 0; i < rows; ++i) {
        // Find the row with maximum element in current column
        int maxRow = i;
        for (int k = i + 1; k < rows; ++k) {
            if (abs(matrix[k][i]) > abs(matrix[maxRow][i])) maxRow = k;
        }
        // Swap the current row with the max row
        swap(matrix[i], matrix[maxRow]);

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < rows; ++k) {
            if (matrix[i][i] == 0) continue;
            double factor = matrix[k][i] / matrix[i][i];
            for (int j = i; j < cols; ++j) {
                matrix[k][j] -= factor * matrix[i][j];
            }
        }
    }
}

// Function to solve the system using back substitution
vector<double> backSubstitution(vector<vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    vector<double> solution(cols - 1, 1.0);

    for (int i = rows - 1; i >= 0; --i) {
        int pivot = -1;
        for (int j = 0; j < cols - 1; ++j) {
            if (matrix[i][j] != 0) {
                pivot = j;
                break;
            }
        }

        if (pivot == -1) continue;

        solution[pivot] = matrix[i][cols - 1];
        for (int j = pivot + 1; j < cols - 1; ++j) {
            solution[pivot] -= matrix[i][j] * solution[j];
        }
        solution[pivot] /= matrix[i][pivot];
    }

    return solution;
}

// Function to find the greatest common divisor
int gcd(int a, int b) {
    return b == 0 ? a : gcd(b, a % b);
}

// Function to convert solution to integer coefficients
vector<int> convertToIntegerCoefficients(const vector<double>& solution) {
    const double EPSILON = 1e-6;
    double minVal = 1.0;
    for (double num : solution) {
        if (abs(num) < EPSILON) continue;
        double frac = abs(num - round(num));
        if (frac > EPSILON) {
            minVal = min(minVal, frac);
        }
    }

    double multiplier = 1.0 / minVal;
    vector<int> intCoeffs;
    for (double num : solution) {
        intCoeffs.push_back(round(num * multiplier));
    }

    // Find the greatest common divisor of all coefficients
    int commonDivisor = abs(intCoeffs[0]);
    for (int i = 1; i < intCoeffs.size(); ++i) {
        commonDivisor = gcd(commonDivisor, abs(intCoeffs[i]));
    }

    // Divide by the GCD to get the simplest form
    for (int& coeff : intCoeffs) {
        coeff /= commonDivisor;
    }

    return intCoeffs;
}

// Function to balance the chemical equation
void balanceEquation(vector<Component>& reactants, vector<Component>& products) {
    auto elements = collectAllElements(reactants, products);
    auto matrix = buildEquationMatrix(elements, reactants, products);
    gaussianElimination(matrix);
    auto solution = backSubstitution(matrix);
    
    solution.push_back(1.0);
    
    auto intCoeffs = convertToIntegerCoefficients(solution);
    // Assign coefficients to reactants and products      
    for (int i = 0; i < reactants.size(); ++i) {
        reactants[i].coefficient = intCoeffs[i];
    }
    for (int i = 0; i < products.size(); ++i) {
        products[i].coefficient = intCoeffs[reactants.size() + i];
    }
}

// Function to print the information about component
void printComponentDetails(const Component& component) {
    cout << "  Component: ";
    if (component.coefficient != 1.0) {
        cout << component.coefficient << " ";
    }
    cout << component.name << endl;
    
    for (const auto& element : component.elements) {
        string elementName = element.first;
        int count = element.second;
        int totalAtoms = count * component.coefficient;
        int oxidationState = component.oxidationStates.at(elementName);
        cout << "    Element: " << elementName 
             << " (Oxidation: " << oxidationState 
             << ", Atoms: " << totalAtoms << ")" << endl;
    }
    cout << "    Total charge: " << component.charge << endl;
}

// Function to get balanced equation as string
string getBalancedEquationString(const vector<Component>& reactants, const vector<Component>& products) {
    stringstream ss;
    // Print reactants
    for (size_t i = 0; i < reactants.size(); ++i) {
        if (i != 0) ss << " + ";
        if (reactants[i].coefficient != 1) ss << reactants[i].coefficient << " ";
        ss << reactants[i].name;
    }
    ss << " = ";
    // Print products
    for (size_t i = 0; i < products.size(); ++i) {
        if (i != 0) ss << " + ";
        if (products[i].coefficient != 1) ss << products[i].coefficient << " ";
        ss << products[i].name;
    }
    return ss.str();
}

void printBalancedEquation(const vector<Component>& reactants, const vector<Component>& products) {
    cout << getBalancedEquationString(reactants, products) << endl;
} 
   
// Print detailed information
void printVerboseOutput(const vector<Component>& reactants, const vector<Component>& products) {
    cout << "Balanced equation: " << getBalancedEquationString(reactants, products) << endl;
    cout << endl;
    
    cout << "Reactants details:" << endl;
    for (const auto& component : reactants) {
        printComponentDetails(component);
    }
    
    cout << "\nProducts details:" << endl;
    for (const auto& component : products) {
        printComponentDetails(component);
    }
}

// Function to normalize the equation (sort components)
string normalizeEquation(const string& equation) {
    size_t equalPos = equation.find('=');
    if (equalPos == string::npos) return equation;
    
    string leftSide = equation.substr(0, equalPos);
    string rightSide = equation.substr(equalPos + 1);
    
    // Separate the components
    auto parseSide = [](const string& side) {
        vector<string> components;
        stringstream ss(side);
        string component;
        while (getline(ss, component, '+')) {
            component.erase(0, component.find_first_not_of(" \t"));
            component.erase(component.find_last_not_of(" \t") + 1);
            if (!component.empty()) {
                components.push_back(component);
            }
        }
        return components;
    };
    
    vector<string> leftComponents = parseSide(leftSide);
    vector<string> rightComponents = parseSide(rightSide);
    
    // Sorting the components
    sort(leftComponents.begin(), leftComponents.end());
    sort(rightComponents.begin(), rightComponents.end());
    
    // Putting it back together
    string normalized;
    for (size_t i = 0; i < leftComponents.size(); ++i) {
        if (i != 0) normalized += " + ";
        normalized += leftComponents[i];
    }
    normalized += " = ";
    for (size_t i = 0; i < rightComponents.size(); ++i) {
        if (i != 0) normalized += " + ";
        normalized += rightComponents[i];
    }
    
    return normalized;
}

// Updated reaction search function
vector<SavedReaction> findReactionsByPartialInput(const string& partialInput) {
    vector<SavedReaction> foundReactions;
    
    string cleanedInput = partialInput;
    cleanedInput.erase(remove(cleanedInput.begin(), cleanedInput.end(), ' '), cleanedInput.end());
    
    // Normalize the input
    string normalizedInput = normalizeEquation(partialInput);
    
    for (const auto& reaction : savedReactions) {
        // Normalizing the stored reaction
        string normalizedSaved = normalizeEquation(reaction.originalEquation);
        
        // We check for a match on either side
        if (normalizedSaved.find(normalizedInput) != string::npos ||
            normalizedInput.find(normalizedSaved) != string::npos) {
            foundReactions.push_back(reaction);
        }
    }
    
    return foundReactions;
}

//Updated reaction save function
void saveReaction(const string& originalEquation, 
                 const vector<Component>& reactants, 
                 const vector<Component>& products) {
    SavedReaction reaction;
    reaction.originalEquation = originalEquation;
    reaction.balancedEquation = getBalancedEquationString(reactants, products);
    reaction.reactants = reactants;
    reaction.products = products;
    
    // We check if there is already such a reaction (normalized)
    string normalizedNew = normalizeEquation(originalEquation);
    for (const auto& saved : savedReactions) {
        if (normalizeEquation(saved.originalEquation) == normalizedNew) {
            return; // reaction already saved
        }
    }
    
    savedReactions.push_back(reaction);
}

// Function to display saved reactions
void displaySavedReactions() {
    if (savedReactions.empty()) {
        cout << "No reactions saved yet." << endl;
        return;
    }
    
    cout << "Saved reactions:" << endl;
    for (size_t i = 0; i < savedReactions.size(); ++i) {
        cout << i + 1 << ". Original: " << savedReactions[i].originalEquation << endl;
        cout << "   Balanced: " << savedReactions[i].balancedEquation << endl;
    }
}

// Function to delete a saved reaction
void deleteSavedReaction() {
    if (savedReactions.empty()) {
        cout << "No reactions to delete." << endl;
        return;
    }
    
    displaySavedReactions();
    cout << "Enter the number of reaction to delete (1-" << savedReactions.size() << "): ";
    
    size_t choice;
    while (!(cin >> choice) || choice < 1 || choice > savedReactions.size()) {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Invalid input. Please enter a number between 1 and " << savedReactions.size() << ": ";
    }
    
    savedReactions.erase(savedReactions.begin() + choice - 1);
    cout << "Reaction deleted successfully!" << endl;
}

// Function to show the main menu
void showMainMenu() {
    cout << "\nChemical Equation Balancer" << endl;
    cout << "1. Balance a new equation" << endl;
    cout << "2. View saved reactions" << endl;
    cout << "3. Delete a saved reaction" << endl;
    cout << "4. Exit" << endl;
    cout << "Enter your choice (1-4): ";
}

void processInteractiveMode() {
    string equation;
    cout << "Enter a chemical equation (e.g., H2 + O2 = H2O) or partial equation (e.g., H2+O2=): ";
    cin.ignore();
    getline(cin, equation);

    size_t equalPos = equation.find('=');
    bool isPartialInput = (equalPos != string::npos) && 
                         (equalPos == equation.length() - 1 || 
                          equation.substr(equalPos + 1).find_first_not_of(" \t") == string::npos);

    if (isPartialInput) {
        auto foundReactions = findReactionsByPartialInput(equation);
        
        if (foundReactions.empty()) {
            cout << "\nNo saved reactions found matching: " << equation << endl;
        } else {
            cout << "\nFound " << foundReactions.size() << " saved reactions:\n";
            for (size_t i = 0; i < foundReactions.size(); ++i) {
                cout << i + 1 << ". " << foundReactions[i].balancedEquation << endl;
            }
            
            cout << "\nEnter number to view details (0 to cancel): ";
            size_t choice;
            while (!(cin >> choice) || choice > foundReactions.size()) {
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
                cout << "Invalid input. Please enter a number between 0 and " << foundReactions.size() << ": ";
            }
            
            if (choice != 0) {
                const auto& selected = foundReactions[choice - 1];
                cout << "\nSelected reaction:\n";
                cout << "Original: " << selected.originalEquation << "\n";
                cout << "Balanced: " << selected.balancedEquation << "\n";
            }
        }
        return;
    }

    if (equalPos == string::npos) {
        cerr << "Invalid equation format." << endl;
        return;
    }

    auto reactants = splitComponents(equation.substr(0, equalPos));
    auto products = splitComponents(equation.substr(equalPos + 1));

    balanceEquation(reactants, products);
    
    cout << "\nBalanced equation:\n";
    printBalancedEquation(reactants, products);
    
    saveReaction(equation, reactants, products);
    cout << "Reaction has been automatically saved.\n";
}

// Function to process the equation in command-line mode
void processInput(int argc, char* argv[]) {
    bool verbose = false;
    string equation;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-v") {
            verbose = true;
        } else {
            // Assume this is the equation (concatenate with spaces if needed)
            if (!equation.empty()) equation += " ";
            equation += arg;
        }
    }

    // If equation wasn't provided as arguments, read from stdin
    if (equation.empty()) {
        cout << "Enter a chemical equation (e.g., H2 + O2 = H2O): ";
        getline(cin, equation);
    }

    size_t equalPos = equation.find('=');
    if (equalPos == string::npos) {
        cerr << "Invalid equation format." << endl;
        return;
    }

    auto reactants = splitComponents(equation.substr(0, equalPos));
    auto products = splitComponents(equation.substr(equalPos + 1));

    balanceEquation(reactants, products);
    
    if (verbose) {
        printVerboseOutput(reactants, products);
    } else {
        printBalancedEquation(reactants, products);
    }

}

int main(int argc, char* argv[]) {
    // Check for --help flag
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--help") {
            displayHelp();
            return 0;
        }
    }

    if (argc > 1) {
        // Command-line mode: join all arguments as the equation
        string equation;
        for (int i = 1; i < argc; ++i) {
            if (i > 1) equation += " ";
            equation += argv[i];
        }
        processInput(argc, argv);
    } else {
        // Interactive mode with menu
        int choice;
        do {
            showMainMenu();
            while (!(cin >> choice) || choice < 1 || choice > 4) {
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
                cout << "Invalid input. Please enter a number between 1 and 4: ";
            }
            
            switch (choice) {
                case 1:
                    processInteractiveMode();
                    break;
                case 2:
                    displaySavedReactions();
                    break;
                case 3:
                    deleteSavedReaction();
                    break;
                case 4:
                    cout << "Exiting program..." << endl;
                    break;
            }
        } while (choice != 4);
    }
    
    return 0;
}