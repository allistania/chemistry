#include "functions.h"

int main(int argc, char* argv[]) {
    loadReactionsFromFile("reactions.txt");

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--help") {
            displayHelp();
            return 0;
        }
    }

    if (argc > 1) {
        processInput(argc, argv);
    } else {
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