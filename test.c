#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <check.h>

#include "lin.h"

int main(int argc, char **argv) {
    struct Matrix *p = createMatrix(1, 2);
    setMatrix(p, 1, 2, (double[1][2]){0.5, 0.5});
    struct Matrix *d = createMatrix(1, 2);
    setMatrix(d, 1, 2, (double[1][2]){1, 0});
    double a = 2;
    double b = 1;

    struct Matrix *s = createMatrix(1, 2);
    struct Matrix *t = createMatrix(1, 2);

    setMatrix(s, 1, 2, (double[1][2]){(p->matrix[0][0])/a, (p->matrix[0][1])/b});
    setMatrix(t, 1, 2, (double[1][2]){(d->matrix[0][0])/a, (d->matrix[0][1])/b});

    // struct Matrix *s2 = mult(s, transpose(s));
    // struct Matrix *t2 = mult(t, transpose(t));

    double s2 = pow(pqNorm(s, 2, 2), 2);
    double t2 = pow(pqNorm(t, 2, 2), 2);

    struct Matrix *st = mult(s, transpose(t));

    printf("s2: %f, t2: %f, st: %f\n", s2, t2, st->matrix[0][0]);
    printf("%f\n", pow(st->matrix[0][0], 2)-(t2*(s2-1)));

    // double lambda = (-st->matrix[0][0] + pow(pow(st->matrix[0][0], 2)-(t2->matrix[0][0]*(s2->matrix[0][0]-1)), 0.5))/(t2->matrix[0][0]);
    double lambda = (double)(-st->matrix[0][0] + pow(pow(st->matrix[0][0], 2)-(t2*(s2-1)), 0.5))/(t2);
    printf("%f\n", lambda);
    return(0);
}
