#define ARRSIZE(A) (sizeof(A) / sizeof(A[0]))

struct Matrix {
    int row;
    int col;
    double **matrix;
};

void printMatrix(struct Matrix *matrix);

void freeMatrix(struct Matrix *A);

void clearMatrix(struct Matrix *A);

struct Matrix *createMatrix(int row, int col);
struct Matrix *createDiagMatrix(int val, int size);

void setMatrix(struct Matrix *A, int row, int col, double values[row][col]);
void equalMatrix(struct Matrix *A, struct Matrix *B);

struct Matrix *mult(struct Matrix *A, struct Matrix *B);
struct Matrix *scalarMult(double scalar, struct Matrix *A);
struct Matrix *add(struct Matrix *A, struct Matrix *B);
struct Matrix *transpose(struct Matrix *A);
double determinant(struct Matrix *A);
struct Matrix *cofactor(struct Matrix *A);
struct Matrix *inverse(struct Matrix *A);

double pqNorm(struct Matrix *A, int p, int q);
double trace(struct Matrix *A);
