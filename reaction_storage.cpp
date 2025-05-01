#include "functions.h"

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

// Load reactions from the file
void loadReactionsFromFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        return;
    }

    savedReactions.clear();

    string line;
    while (getline(file, line)) {
        vector<string> parts;
        stringstream ss(line);
        string part;
        while (getline(ss, part, '|')) {
            part.erase(0, part.find_first_not_of(" \t"));
            part.erase(part.find_last_not_of(" \t") + 1);
            parts.push_back(part);
        }

        if (parts.size() != 4) continue;

        string originalEquation = parts[0];
        string balancedEquation = parts[1];
        string reactantCoeffsStr = parts[2];
        string productCoeffsStr = parts[3];

        vector<double> reactantCoeffs;
        stringstream rcStream(reactantCoeffsStr);
        string rc;
        while (getline(rcStream, rc, ',')) {
            reactantCoeffs.push_back(stod(rc));
        }

        vector<double> productCoeffs;
        stringstream pcStream(productCoeffsStr);
        string pc;
        while (getline(pcStream, pc, ',')) {
            productCoeffs.push_back(stod(pc));
        }

        size_t equalPos = originalEquation.find('=');
        if (equalPos == string::npos) continue;

        vector<Component> reactants = splitComponents(originalEquation.substr(0, equalPos));
        vector<Component> products = splitComponents(originalEquation.substr(equalPos + 1));

        if (reactants.size() != reactantCoeffs.size() || products.size() != productCoeffs.size()) {
            continue;
        }

        for (size_t i = 0; i < reactants.size(); ++i) {
            reactants[i].coefficient = reactantCoeffs[i];
        }
        for (size_t i = 0; i < products.size(); ++i) {
            products[i].coefficient = productCoeffs[i];
        }

        SavedReaction reaction;
        reaction.originalEquation = originalEquation;
        reaction.balancedEquation = balancedEquation;
        reaction.reactants = reactants;
        reaction.products = products;

        savedReactions.push_back(reaction);
    }

    file.close();
}

// Save reactions to the file
void saveReactionsToFile(const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error saving reactions to file." << endl;
        return;
    }

    for (const auto& reaction : savedReactions) {
        file << reaction.originalEquation << " | " << reaction.balancedEquation;

        file << " | ";
        for (size_t i = 0; i < reaction.reactants.size(); ++i) {
            if (i > 0) file << ",";
            file << reaction.reactants[i].coefficient;
        }

        file << " | ";
        for (size_t i = 0; i < reaction.products.size(); ++i) {
            if (i > 0) file << ",";
            file << reaction.products[i].coefficient;
        }

        file << endl;
    }

    file.close();
}

// Updated reaction save function
void saveReaction(const string& originalEquation, 
                 const vector<Component>& reactants, 
                 const vector<Component>& products) {
    SavedReaction reaction;
    reaction.originalEquation = originalEquation;
    reaction.balancedEquation = getBalancedEquationString(reactants, products);
    reaction.reactants = reactants;
    reaction.products = products;

    string normalizedNew = normalizeEquation(originalEquation);
    bool exists = false;
    for (const auto& saved : savedReactions) {
        if (normalizeEquation(saved.originalEquation) == normalizedNew) {
            exists = true;
            break;
        }
    }

    if (!exists) {
        savedReactions.push_back(reaction);
        saveReactionsToFile("reactions.txt");
    }
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
    saveReactionsToFile("reactions.txt");
    cout << "Reaction deleted successfully!" << endl;
}
