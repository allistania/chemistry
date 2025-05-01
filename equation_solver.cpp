#include "functions.h"
#include <numeric>

using namespace std;

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
    vector<double> solution(cols - 1, 0.0);

    // Let's assume that the last variable is free. (x_n = 1)
    solution[cols - 2] = 1.0;

    for (int i = rows - 1; i >= 0; --i) {
        double sum = 0.0;
        int pivot = -1;

        for (int j = 0; j < cols - 1; ++j) {
            if (abs(matrix[i][j]) > 1e-12) {
                if (pivot == -1) pivot = j;
                else sum += matrix[i][j] * solution[j];
            }
        }

        if (pivot != -1)
            solution[pivot] = (matrix[i][cols - 1] - sum) / matrix[i][pivot];
    }

    return solution;
}
// Function to find the greatest common divisor
int gcd(int a, int b) {
    return b == 0 ? a : gcd(b, a % b);
}

// Fraction structure
struct Fraction {
    int numerator;
    int denominator;

    Fraction(double value, double epsilon = 1e-8) {
        int sign = value < 0 ? -1 : 1;
        value = fabs(value);

        int bestDen = 1;
        double minError = fabs(value - round(value));

        for (int d = 1; d <= 1000; ++d) {
            int n = round(value * d);
            double error = fabs(value - double(n) / d);
            if (error < epsilon) {
                numerator = sign * n;
                denominator = d;
                return;
            }
            if (error < minError) {
                bestDen = d;
                minError = error;
            }
        }

        numerator = sign * round(value * bestDen);
        denominator = bestDen;
    }
};

// Function to convert fractional coefficients to integer coefficients
vector<int> convertToIntegerCoefficients(const vector<double>& solution) {
    vector<Fraction> fracs;
    vector<int> denominators;

    for (double val : solution) {
        Fraction f(val);
        fracs.push_back(f);
        denominators.push_back(f.denominator);
    }

    // Find the least common multiple of the denominators
    int commonDen = denominators[0];
    for (int i = 1; i < denominators.size(); ++i) {
        commonDen = lcm(commonDen, denominators[i]);
    }

    // We get integer coefficients
    vector<int> intCoeffs;
    for (const auto& frac : fracs) {
        intCoeffs.push_back(frac.numerator * (commonDen / frac.denominator));
    }

    // We reduce by GCD
    int commonGCD = abs(intCoeffs[0]);
    for (int i = 1; i < intCoeffs.size(); ++i) {
        commonGCD = gcd(commonGCD, abs(intCoeffs[i]));
    }

    for (int& val : intCoeffs) {
        val /= commonGCD;
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
