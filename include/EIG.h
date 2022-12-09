#pragma once

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <random>

#include "Matrix.h"
#include "Vector.h"
#include "QR.h"

// Define error codes.
constexpr int EIG_MATRIXNOTSQUARE = -1;
constexpr int EIG_MAXITERATIONSEXCEEDED = -2;
constexpr int EIG_MATRIXNOTSYMMETRIC = -3;

// Function to estimate (real) eigenvalues using QR decomposition.
/* Note that this is only valid for matrices that have ALL real
	eigenvalues. The only matrices that are guaranteed to have only
	real eigenvalues are symmetric matrices. Therefore, this function
	is only guaranteed to work with symmetric matrices. */
template <typename T>
int EigQR(const Matrix2<T> &inputMatrix, std::vector<T> &eigenValues)
{
    // Make a copy of the input matrix.
    Matrix2<T> A = inputMatrix;

    // Verify that the input matrix is square.
    if (!A.IsSquare())
        return EIG_MATRIXNOTSQUARE;

    // Verify that the matrix is symmetric.
    if (!A.IsSymmetric())
        return EIG_MATRIXNOTSYMMETRIC;

    // The number of eigenvalues is equal to the number of rows.
    int numRows = A.GetNumRows();

    // Create an identity matrix of the same dimensions.
    Matrix2<T> identityMatrix(numRows, numRows);
    identityMatrix.SetToIdentity();

    // Create matrices to store Q and R.
    Matrix2<T> Q (numRows, numRows);
    Matrix2<T> R (numRows, numRows);

    // Loop through each iteration.
    int maxIterations = 10e3;
    int iterationCount = 0;
    bool continueFlag = true;
    while ((iterationCount < maxIterations) && continueFlag)
    {
        // Compute the QR decomposition of A.
        int returnValue = QR<T>(A, Q, R);

        // Compute the next value of A as the product of R and Q.
        A = R * Q;

        /* Check if A is now close enough to being upper-triangular.
			We can do this using the IsRowEchelon() function from the
			Matrix2 class. */
        if (A.IsRowEchelon())
            continueFlag = false;

        // Increment iterationCount.
        iterationCount++;
    }

    // At this point, the eigenvalues should be the diagonal elements of A.
    for (int i=0; i<numRows; ++i)
        eigenValues.push_back(A.GetElement(i,i));

    // Set the return status accordingly.
    if (iterationCount == maxIterations)
        return EIG_MAXITERATIONSEXCEEDED;
    else
        return 0;

}

// Function to perform inverse power iteration method.
template <typename T>
int InvPIt(const Matrix2<T> &inputMatrix, const T &eigenValue, Vector<T> &eigenVector)
{
    // Make a copy of the input matrix.
    Matrix2<T> A = inputMatrix;

    // Verify that the input matrix is square.
    if (!A.IsSquare())
        return EIG_MATRIXNOTSQUARE;

    // Setup a random number generator.
    std::random_device myRandomDevice;
    std::mt19937 myRandomGenerator(myRandomDevice());
    std::uniform_int_distribution<int> myDistribution(1.0, 10.0);

    /* The number of eigenvectors and eigenvalues that we will compute will be
		equal to the number of rows in the input matrix. */
    int numRows = A.GetNumRows();

    // Create an identity matrix of the same dimensions.
    Matrix2<T> identityMatrix(numRows, numRows);
    identityMatrix.SetToIdentity();

    // Create an initial vector, v.
    Vector<T> v(numRows);
    for (int i=0; i<numRows; ++i)
        v.SetElement(i, static_cast<T>(myDistribution(myRandomGenerator)));

    // Iterate.
    int maxIterations = 100;
    int iterationCount = 0;
    T deltaThreshold = static_cast<T>(1e-9);
    T delta = static_cast<T>(1e6);
    Vector<T> prevVector(numRows);
    Matrix2<T> tempMatrix(numRows, numRows);

    while ((iterationCount < maxIterations) && (delta > deltaThreshold))
    {
        // Store a copy of the current working vector to use for computing delta.
        prevVector = v;

        // Compute the next value of v.
        tempMatrix = A - (eigenValue * identityMatrix);
        tempMatrix.Inverse();
        v = tempMatrix * v;
        v.Normalize();

        // Compute delta.
        delta = (v - prevVector).norm();

        // Increment iteration count.
        iterationCount++;
    }

    // Return the estimated eigenvector.
    eigenVector = v;

    // Set the return status accordingly.
    if (iterationCount == maxIterations)
        return EIG_MAXITERATIONSEXCEEDED;
    else
        return 0;

}
