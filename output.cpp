#include "functions.h"

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

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-v") {
            verbose = true;
        } else {
            if (!equation.empty()) equation += " ";
            equation += arg;
        }
    }

    if (equation.empty()) {
        cout << "Enter a chemical equation (e.g., H2 + O2 = H2O): ";
        getline(cin, equation);
    }

    bool isPartial = false;
    size_t equalPos = equation.find('=');
    if (equalPos != string::npos) {
        string rightSide = equation.substr(equalPos + 1);
        string leftSide = equation.substr(0, equalPos);
        if (rightSide.find_first_not_of(" \t") == string::npos ||
            leftSide.find_first_not_of(" \t") == string::npos) {
            isPartial = true;
        }
    }

    if (isPartial) {
        vector<SavedReaction> found = findReactionsByPartialInput(equation);
        if (found.empty()) {
            cout << "No matching reactions found." << endl;
        } else {
            cout << "Found " << found.size() << " reactions:" << endl;
            for (const auto& r : found) {
                cout << r.balancedEquation << endl;
            }
        }
        return;
    }

    equalPos = equation.find('=');
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

    saveReaction(equation, reactants, products);
}
