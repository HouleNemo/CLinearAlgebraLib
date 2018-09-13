#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lin.h"

/**
 Prints matrix out to stdout.
 */
void printMatrix(struct Matrix *matrix) {
    for (int i = 0; i < matrix->row; i++) {
        for (int j = 0; j < matrix->col; j++) {
            fprintf(stdout, "%f ", matrix->matrix[i][j]);
        }
        fprintf(stdout, "\n");
    }
}

/**
 Free matrix A.
 */
void freeMatrix(struct Matrix *A) {
    for (int i = 0; i < A->row; i++) {
        free(A->matrix[i]);
    }
    free(A->matrix);
    free(A);
}

/**
 Set all elements in A to 0.
 */
void clearMatrix(struct Matrix *A) {
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            A->matrix[i][j] = 0;
        }
    }
}

/**
 Create matrix with given # of rows and # of columns, all elements set to 0.
 */
struct Matrix *createMatrix(int row, int col) {
    struct Matrix * R = (struct Matrix *) malloc(sizeof(struct Matrix));
    R->row = row;
    R->col = col;
    R->matrix = (double **) malloc(row*sizeof(double *));
    for (int i = 0; i < row; i++) {
        R->matrix[i] = (double *) malloc(col*sizeof(double));
    }
    clearMatrix(R);
    return R;
}
/**
 Create square matrix with # of rows and # of columns equal to size. All
 elements on diagonal are set to val, everything else is 0.
 */
struct Matrix *createDiagMatrix(int val, int size) {
    struct Matrix *R = createMatrix(size, size);
    for (int i = 0; i < size; i++) {
        R->matrix[i][i] = val;
    }
    return R;
}

/**
 Store 2D array of values into A, row and col must match with # of rows and #
 of columns in A
 */
void setMatrix(struct Matrix *A, int row, int col, double values[row][col]) {
    if (A->row != row || A->col != col) {
        fprintf(stderr, "Invalid matrix and values size in setMatrix\n");
        return;
    }

    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            A->matrix[i][j] = values[i][j];
        }
    }
}
/**
 Deep copy elements of B into A, both must be of the same size
 */
void equalMatrix(struct Matrix *A, struct Matrix *B) {
    if (A->row != B->row || A->col != B->col) {
        fprintf(stderr, "Invalid matrix size for equalMatrix\n");
        return;
    }

    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            A->matrix[i][j] = B->matrix[i][j];
        }
    }
}

/**
 Returns matrix multiplication of AB, # of columns in A must match # of rows
 in B
 */
struct Matrix *mult(struct Matrix *A, struct Matrix *B) {

    if (A->col != B->row) {
        fprintf(stderr, "Invalid matrix sizes for multiplying.\n");
        return NULL;
    }

    struct Matrix *R = createMatrix(A->row, B->col);
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            for (int k = 0; k < B->col; k++) {
                R->matrix[i][k] += A->matrix[i][j] * B->matrix[j][k];
            }
        }
    }
    return R;
}
/**
 Returns scalar multiplication of scalar*A
 */
struct Matrix *scalarMult(double scalar, struct Matrix *A) {
    struct Matrix *R = createMatrix(A->row, A->col);
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            R->matrix[i][j] = scalar * A->matrix[i][j];
        }
    }
    return R;
}
/**
 Returns matrix addition of A+B, both must be of the same size
 */
struct Matrix *add(struct Matrix *A, struct Matrix *B) {

    if (A->row != B->row|| A->col != B->col) {
      fprintf(stderr, "Invalid matrix sizes for adding.\n");
      return NULL;
    }

    struct Matrix *R = createMatrix(A->row, A->col);
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            R->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
    }
    return R;
}
/**
 Returns the transpose of A
 */
struct Matrix *transpose(struct Matrix *A) {
    struct Matrix *R = createMatrix(A->col, A->row);
    for (int i = 0; i < A->row; i++) {
        for (int j = 0; j < A->col; j++) {
            R->matrix[j][i] = A->matrix[i][j];
        }
    }
    return R;
}
/**
 Calculates the determinant of A, must be a square matrix
 */
double determinant(struct Matrix *A) {
    /* Adapted from:
    https://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html */
    if (A->row != A->col) {
        fprintf(stderr, "Matrix is not a square matrix, cannot find determinant");
        return 0;
    }

    double det = 1;
    int n = A->row;
    struct Matrix *subA = createMatrix(n-1, n-1);

    if (n == 1) {
        det = A->matrix[0][0];
    } else if (n == 2) {
        det = A->matrix[0][0]*A->matrix[1][1] - A->matrix[0][1]*A->matrix[1][0];
    } else {
        for (int j1=0;j1<n;j1++) {
            clearMatrix(subA);
            for (int i=1;i<n;i++) {
                int j2 = 0;
                for (int j=0;j<n;j++) {
                    if (j == j1){continue;};
                    subA->matrix[i-1][j2] = A->matrix[i][j];
                    j2++;
                }
            }
            det += pow(-1.0,j1+2.0) * A->matrix[0][j1] * determinant(subA);
      }
    }
    freeMatrix(subA);

    return det;
}
/**
 Calculates the CoFactor matrix of A
 */
struct Matrix *cofactor(struct Matrix *A) {
    /* Adapted from:
    https://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html */
    int i,j,ii,jj,i1,j1;
    double det;
    int n = A->row;
    struct Matrix *temp = createMatrix(n-1, n-1);
    struct Matrix *R = createMatrix(n, n);

    if (n != A->col) {
        fprintf(stderr, "Invalid matrix sizes for cofactor\n");
        return NULL;
    }

    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) {

            /* Form the adjoint a_ij */
            i1 = 0;
            for (ii=0;ii<n;ii++) {
                if (ii == i){continue;};
                j1 = 0;
                for (jj=0;jj<n;jj++) {
                    if (jj == j){continue;};
                    temp->matrix[i1][j1] = A->matrix[ii][jj];
                    j1++;
                }
                i1++;
            }

            /* Calculate the determinant */
            det = determinant(temp);

            /* Fill in the elements of the cofactor */
            R->matrix[i][j] = pow(-1.0,i+j+2.0) * det;
        }
    }

    freeMatrix(temp);
    return R;
}
/**
 Calculates the inverse of A, must be a square matrix
 */
struct Matrix *inverse(struct Matrix *A) {
    if (A->row != A->col) {
        fprintf(stderr, "Invalid matrix sizes for inverse\n");
        return NULL;
    }

    double det = determinant(A);
    struct Matrix *cf = cofactor(A);
    struct Matrix *adjoint = transpose(cf);
    struct Matrix *R = scalarMult(1/det, adjoint);

    freeMatrix(cf);
    freeMatrix(adjoint);
    return R;
}

/**
 Calculates the "entrywise" matrix norm
 https://en.wikipedia.org/wiki/Matrix_norm#L2,1_and_Lp,q_norms
 */
double pqNorm(struct Matrix *A, int p, int q) {
    double sum = 0;
    for (int j = 0; j < A->col; j++) {
        double isum = 0;
        for (int i = 0; i < A->row; i++) {
            isum += pow(fabs(A->matrix[i][j]), p);
        }
        sum += pow(isum, (double)q/p);
    }
    return pow(sum, (double)1/q);
}

/**
 Calculates the trace of a square matrix
 https://en.wikipedia.org/wiki/Trace_(linear_algebra)
 */
double trace(struct Matrix *A) {
    if (A->row != A->col) {
        fprintf(stderr, "Invalid matrix sizes for trace\n");
        return 0;
    }

    double sum = 0;
    for (int i = 0; i < A->row; i++) {
        sum += A->matrix[i][i];
    }
    return sum;
}
