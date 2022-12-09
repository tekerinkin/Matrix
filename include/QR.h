#pragma once

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "Matrix.h"
#include "Vector.h"


constexpr int QR_MATRIXNOTSQUARE = -1;


template <typename T>
int QR(const Matrix2<T> &A, Matrix2<T> &Q, Matrix2<T> &R)
{

    // Make a copy of the input matrix.
    Matrix2<T> inputMatrix = A;

    // Verify that the input matrix is square.
    if (!inputMatrix.IsSquare())
        return QR_MATRIXNOTSQUARE;

    // Determine the number of columns (and rows, since the matrix is square).
    int numCols = inputMatrix.GetNumCols();

    // Create a vector to store the P matrices for each column.
    std::vector<Matrix2<T>> Plist;

    // Loop through each column.
    for (int j=0; j<(numCols-1); ++j)
    {
        // Create the a1 and b1 vectors.
        // a1 is the column vector from A.
        // b1 is the vector onto which we wish to reflect a1.
        Vector<T> a1 (numCols-j);
        Vector<T> b1 (numCols-j);
        for (int i=j; i<numCols; ++i)
        {
            a1.SetElement(i-j, inputMatrix.GetElement(i,j));
            b1.SetElement(i-j, static_cast<T>(0.0));
        }
        b1.SetElement(0, static_cast<T>(1.0));

        // Compute the norm of the a1 vector.
        T a1norm = a1.norm();

        // Compute the sign we will use.
        int sgn = -1;
        if (a1.GetElement(0) < static_cast<T>(0.0))
            sgn = 1;

        // Compute the u-vector.
        Vector<T> u = a1 - (sgn * a1norm * b1);

        // Compute the n-vector.
        Vector<T> n = u.Normalized();

        // Convert n to a matrix so that we can transpose it.
        Matrix2<T> nMat (numCols-j, 1);
        for (int i=0; i<(numCols-j); ++i)
            nMat.SetElement(i, 0, n.GetElement(i));

        // Transpose nMat.
        Matrix2<T> nMatT = nMat.Transpose();

        // Create an identity matrix of the appropriate size.
        Matrix2<T> I (numCols-j, numCols-j);
        I.SetToIdentity();

        // Compute Ptemp.
        Matrix2<T> Ptemp = I - static_cast<T>(2.0) * nMat * nMatT;

        // Form the P matrix with the original dimensions.
        Matrix2<T> P (numCols, numCols);
        P.SetToIdentity();
        for (int row=j; row<numCols; ++row)
        {
            for (int col=j; col<numCols; ++col)
            {
                P.SetElement(row, col, Ptemp.GetElement(row-j, col-j));
            }
        }

        // Store the result into the Plist vector.
        Plist.push_back(P);

        // Apply this transform matrix to inputMatrix and use this result
        // next time through the loop.
        inputMatrix = P * inputMatrix;
    }

    // Compute Q.
    Matrix2<T> Qmat = Plist.at(0);
    for (int i=1; i<(numCols-1); ++i)
    {
        Qmat = Qmat * Plist.at(i).Transpose();
    }

    // Return the Q matrix.
    Q = Qmat;

    // Compute R.
    int numElements = Plist.size();
    Matrix2<T> Rmat = Plist.at(numElements-1);
    for (int i=(numElements-2); i>=0; --i)
    {
        Rmat = Rmat * Plist.at(i);
    }
    Rmat = Rmat * A;

    // And return the R matrix.
    R = Rmat;

}
