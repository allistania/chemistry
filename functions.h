#pragma once

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
#include <fstream>

using namespace std;

// Data structures
struct Component {
    string name;
    map<string, int> elements;
    map<string, int> oxidationStates;
    int charge;
    double coefficient;
    bool isOrganic;
};

struct SavedReaction {
    string originalEquation;
    string balancedEquation;
    vector<Component> reactants;
    vector<Component> products;
};

// Global variables
extern vector<SavedReaction> savedReactions;
extern unordered_map<string, int> fixedOxidationStates;
extern unordered_map<string, vector<int>> variableOxidationStates;

// Function prototypes
void displayHelp();
vector<pair<string, int>> parseElementCounts(const string& componentStr);
Component parseComponent(const string& componentStr);
vector<Component> splitComponents(const string& side);
vector<string> collectAllElements(const vector<Component>& reactants, const vector<Component>& products);
vector<vector<double>> buildEquationMatrix(const vector<string>& elements, const vector<Component>& reactants, const vector<Component>& products);
void gaussianElimination(vector<vector<double>>& matrix);
vector<double> backSubstitution(vector<vector<double>>& matrix);
vector<int> convertToIntegerCoefficients(const vector<double>& solution);
void balanceEquation(vector<Component>& reactants, vector<Component>& products);
void printComponentDetails(const Component& component);
string getBalancedEquationString(const vector<Component>& reactants, const vector<Component>& products);
void printBalancedEquation(const vector<Component>& reactants, const vector<Component>& products);
void printVerboseOutput(const vector<Component>& reactants, const vector<Component>& products);
string normalizeEquation(const string& equation);
vector<SavedReaction> findReactionsByPartialInput(const string& partialInput);
void saveReaction(const string& originalEquation, const vector<Component>& reactants, const vector<Component>& products);
void displaySavedReactions();
void deleteSavedReaction();
void showMainMenu();
void processInteractiveMode();
void processInput(int argc, char* argv[]);
void loadReactionsFromFile(const string& filename);
void saveReactionsToFile(const string& filename);