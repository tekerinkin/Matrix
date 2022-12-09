#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <fstream>

#include "../include/Matrix.h"
#include "../include/Vector.h"
#include "../include/EIG.h"

using namespace std;

// A simple function to print a matrix to stdout.
template <class T>
void PrintMatrix(Matrix2<T> matrix)
{
    int nRows = matrix.GetNumRows();
    int nCols = matrix.GetNumCols();
    for (int row = 0; row<nRows; ++row)
    {
        for (int col = 0; col<nCols; ++col)
        {
            cout << std::fixed << std::setprecision(3) << matrix.GetElement(row, col) << "  ";
        }
        cout << endl;
    }
}

// A simple function to print a vector to stdout.
template <class T>
void PrintVector(Vector<T> inputVector)
{
    int nRows = inputVector.GetNumDims();
    for (int row = 0; row<nRows; ++row)
    {
        cout << std::fixed << std::setprecision(6) << inputVector.GetElement(row) << endl;
    }
}



int main() {
    cout << "**********************************************" << endl;
    cout << "Testing eigenvalue and eigenvector code." << endl;
    cout << "Eigenvalues by QR decomposition." << endl;
    cout << "**********************************************" << endl;
    cout << endl;

    {
        cout << "Testing with a simple (non-symmetric) 3x3 matrix:" << endl;
        std::vector<double> simpleData = {0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
        Matrix2<double> testMatrix(3, 3, simpleData);
        PrintMatrix(testMatrix);
        cout << endl;

        // Compute the eigenvalues.
        std::vector<double> eigenValues;
        int returnStatus = EigQR(testMatrix, eigenValues);

        if (returnStatus == EIG_MAXITERATIONSEXCEEDED)
            cout << ">>> Maximum iterations exceeded <<<" << endl;

        if (returnStatus == EIG_MATRIXNOTSYMMETRIC)
            cout << ">>> Matrix not symmetric. <<<" << endl;

        // Display the eigenvalues.
        cout << "The estimated eigenvalues are:" << endl;
        for (auto currentValue : eigenValues)
            cout << std::setprecision(6) << currentValue << " ";

        cout << endl << endl;
    }

    {
        cout << "Testing with a simple (symmetric) 3x3 matrix:" << endl;
        std::vector<double> simpleData = {6.0, 5.5, -1.0, 5.5, 1.0, -2.0, -1.0, -2.0, -3.0};
        Matrix2<double> testMatrix(3, 3, simpleData);
        PrintMatrix(testMatrix);
        cout << endl;

        // Compute the eigenvalues.
        std::vector<double> eigenValues;
        int returnStatus = EigQR(testMatrix, eigenValues);

        if (returnStatus == EIG_MAXITERATIONSEXCEEDED)
            cout << ">>> Maximum iterations exceeded <<<" << endl;

        if (returnStatus == EIG_MATRIXNOTSYMMETRIC)
            cout << ">>> Matrix not symmetric. <<<" << endl;

        // Display the eigenvalues.
        cout << "The estimated eigenvalues are:" << endl;
        for (auto currentValue : eigenValues)
            cout << std::setprecision(6) << currentValue << " ";

        cout << endl << endl;

        // Setup a vector for the eigenvector.
        Vector<double> eigenVector(3);

        // Loop through each eigenvalue.
        for (auto currentValue : eigenValues)
        {
            cout << "Estimated eigenvector for eigenvalue " << currentValue << " = " << endl;
            int returnStatus = InvPIt<double>(testMatrix, currentValue, eigenVector);
            PrintVector(eigenVector);

            if (returnStatus == EIG_MAXITERATIONSEXCEEDED)
                cout << "*** Maximum iterations exceeded ***" << endl;

            cout << endl;
        }

        cout << endl << endl;
    }

    {
        cout << "Testing with an example that should have complex eigenvalues:" << endl;
        std::vector<double> simpleData = {4.0, -6.0, 8.0, 7.0, 9.0, -5.0, 9.0, -6.0, -4.0};
        Matrix2<double> testMatrix(3, 3, simpleData);
        PrintMatrix(testMatrix);
        cout << endl;

        // Compute the eigenvalues.
        std::vector<double> eigenValues;
        int returnStatus = EigQR(testMatrix, eigenValues);

        if (returnStatus == EIG_MATRIXNOTSYMMETRIC)
            cout << ">>> Matrix not symmetric. <<<" << endl;

        if (returnStatus == EIG_MAXITERATIONSEXCEEDED)
            cout << ">>> Maximum iterations exceeded <<<" << endl;

        // Display the eigenvalues.
        cout << "The estimated eigenvalues are:" << endl;
        for (auto currentValue : eigenValues)
            cout << std::setprecision(6) << currentValue << " ";

        cout << endl << endl;

        // Setup a vector for the eigenvector.
        Vector<double> eigenVector(3);

        // Loop through each eigenvalue.
        for (auto currentValue : eigenValues)
        {
            cout << "Estimated eigenvector for eigenvalue " << currentValue << " = " << endl;
            int returnStatus = InvPIt<double>(testMatrix, currentValue, eigenVector);
            PrintVector(eigenVector);

            if (returnStatus == EIG_MAXITERATIONSEXCEEDED)
                cout << "*** Maximum iterations exceeded ***" << endl;

            cout << endl;
        }

        cout << endl << endl;
    }

    {
        cout << "Testing with random matrix: " << endl;
        cout << "Input size: " << endl;
        int n;
        cin >> n;

        cout << "Matrix size: " << n << "x" << n << endl;
        std::vector<double> simpleData(n*n);

        for(size_t i = 0; i < n; ++i) {
            for(size_t j = 0; j < n; ++j) {
                simpleData[i*n + j] = rand()%100 + 10;
                simpleData[j*n + i] = simpleData[i*n + j];
            }
        }

        Matrix2<double> testMatrix(n, n, simpleData);
        PrintMatrix(testMatrix);
        cout << endl;

        // Compute the eigenvalues.
        std::vector<double> eigenValues;
        int returnStatus = EigQR(testMatrix, eigenValues);

        if (returnStatus == EIG_MATRIXNOTSYMMETRIC)
            cout << ">>> Matrix not symmetric. <<<" << endl;

        if (returnStatus == EIG_MAXITERATIONSEXCEEDED)
            cout << ">>> Maximum iterations exceeded <<<" << endl;

        // Display the eigenvalues.
        cout << "The estimated eigenvalues are:" << endl;
        for (auto currentValue : eigenValues)
            cout << std::setprecision(6) << currentValue << " ";

        cout << endl << endl;

        // Setup a vector for the eigenvector.
        Vector<double> eigenVector(3);

        // Loop through each eigenvalue.
        for (auto currentValue : eigenValues)
        {
            cout << "Estimated eigenvector for eigenvalue " << currentValue << " = " << endl;
            int returnStatus = InvPIt<double>(testMatrix, currentValue, eigenVector);
            PrintVector(eigenVector);

            if (returnStatus == EIG_MAXITERATIONSEXCEEDED)
                cout << "*** Maximum iterations exceeded ***" << endl;

            cout << endl;
        }

        cout << endl << endl;
    }

    return 0;
}