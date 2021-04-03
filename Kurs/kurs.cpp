#include <stdio.h>
#include "quanc8.h"
#include "matrix.h"
#include "rkf45.h"
#include <math.h>
#include <iostream>

using namespace std;

// Подынтергальная функция
double function(double t)
{
    return 1.0 / (1.0 + cos(t));
}

//Решение интеграла A
double integral() 
{
    const double PI = acos(-1.0);
    double a = 0;             // нижний предел интегрирования
    double b = PI / 2.0;      // верхний предел интегрирования
    double abserr = 1.0e-14;  // абсолютная погрешность
    double relerr = 1.0e-14;  // относительная погрешность

    double errest;            // оценка погрешности
    int nofun;                // количество вычислений подынтегральной функции
    double flag;              // индикатор надежности

    // считаем при помощи QUANC8
    double f;
    quanc8(function, a, b, abserr, relerr, &f, &errest, &nofun, &flag);

    // вывод
    //printf("QUANC8 = ");
    //printf("%6.4f ", f);
    //cout << endl;
    return f;
}

void printMatrix(MATRIX(d), int n)
{
    for (int i = 0; i < n; i++)
    {
        if (i == 0) printf("  %c", 218);
        else
            if (i == (n - 1)) printf("  %c", 192);
            else
                printf("  %c", 179);
        for (int j = 0; j < n; j++)
        {
            printf("%10.4f", d[i][j]);
        }
        if (i == 0) printf("   %c", 191);
        else
            if (i == (n - 1)) printf("   %c", 217);
            else
                printf("   %c", 179);
        printf("\n");
    }
}

void printVector(double* vector, int n)
{
    for (int i = 0; i < n; i++)
    {
        if (i == 0) printf("  %c", 218);
        else
            if (i == (n - 1)) printf("  %c", 192);
            else
                printf("  %c", 179);
        printf("%10.4f", vector[i]);
        if (i == 0) printf("   %c", 191);
        else
            if (i == (n - 1)) printf("   %c", 217);
            else
                printf("   %c", 179);
        printf("\n");
    }
}

double* solveSystemOfLinearEquations(MATRIX(input), double* vector, double* cond, int size)
{
    double* work = new double[size];
    int* ipvt = new int[size];
    decomp(size, input, cond, ipvt, work);
    solve(size, input, vector, ipvt);
    double* result = vector;
    return result;
}

void MyFunc(float t, float* y, float* dy)
{
    // y" + x * y' - y / (2 * x) = A
    // y" = - x * y' + y / (2 * x) + A
    dy[0] = y[1];
    double A = integral();
    dy[1] = -t * y[1] + y[0] / (2 * t) + A;
}

void main(void)
{
    //Решение методом конечных разностей при h = 0.1

    cout << "  h = 0.1" << endl << endl;
    MATRIX(matrix);
    double cond;
    double vector[4];

    matrix[0][0] = 2.9;
    matrix[0][1] = -4.0;
    matrix[0][2] = 1.0;
    matrix[0][3] = 0.0;

    cout << "   2.9 * y_0 -   4.0    * y_1 +   1.0    * y_2 +   0.0 * y_3 = -0.1" << endl;

    matrix[1][0] = 89.5;
    matrix[1][1] = -200.2381;
    matrix[1][2] = 110.5;
    matrix[1][3] = 0.0;

    cout << "  89.5 * y_0 - 200.2381 * y_1 + 110.5    * y_2 +   0.0 * y_3 =  1.0" << endl;

    matrix[2][0] = 0.0;
    matrix[2][1] = 89.0;
    matrix[2][2] = -200.2273;
    matrix[2][3] = 111.0;

    cout << "   0.0 * y_0 +  89.0    * y_1 - 200.2273 * y_2 + 111.0 * y_3 =  1.0" << endl;

    matrix[3][0] = 0.0;
    matrix[3][1] = 0.0;
    matrix[3][2] = 0.0;
    matrix[3][3] = 1.0;

    cout << "   0.0 * y_0 +   0.0    * y_1 +   0.0    * y_2 +   1.0 * y_3 =  2.1599" << endl << endl;

    printMatrix(matrix, 4);
    cout << endl;

    vector[0] = -0.1;
    vector[1] = integral();
    vector[2] = integral();
    vector[3] = 2.1599;

    printVector(vector, 4);
    cout << endl;

    double* y0_3 = solveSystemOfLinearEquations(matrix, vector, &cond, 4);

    cout << "  cond = ";
    printf("%9.4f", cond);
    cout << endl << endl;

    cout << "  ";
    for (int i = 0; i < 4; i++) {
        cout << "x_" << i << " = ";
        printf("%3.1f      ", 2.0 + i * 0.1);
    }
    cout << endl;

    cout << "  ";
    for (int i = 0; i < 4; i++) {
        cout << "y_" << i << " = ";
        printf("%6.4f   ", y0_3[i]);
    }
    cout << endl << endl;

    //Решение методом конечных разностей при h = 0.05

    cout << "  h = 0.05" << endl << endl;
    MATRIX(matrix2);
    double vector2[7];

    matrix2[0][0] = 2.95;
    matrix2[0][1] = -4.0;
    matrix2[0][2] = 1.0;
    matrix2[0][3] = 0.0;
    matrix2[0][4] = 0.0;
    matrix2[0][5] = 0.0;
    matrix2[0][6] = 0.0;

    cout << "  2.95 * y_0 - 4.0 * y_1 + 1.0 * y_2 = -0.05" << endl;

    matrix2[1][0] = 379.5;
    matrix2[1][1] = -800.2439;
    matrix2[1][2] = 420.5;
    matrix2[1][3] = 0.0;
    matrix2[1][4] = 0.0;
    matrix2[1][5] = 0.0;
    matrix2[1][6] = 0.0;

    cout << "  379.5 * y_0 - 800.2439 * y_1 + 420.5 * y_2 = 1.0" << endl;

    matrix2[2][0] = 0.0;
    matrix2[2][1] = 379.0;
    matrix2[2][2] = -800.2381;
    matrix2[2][3] = 421.0;
    matrix2[2][4] = 0.0;
    matrix2[2][5] = 0.0;
    matrix2[2][6] = 0.0;

    cout << "  379.0 * y_1 - 800.2381 * y_2 + 421.0 * y_3 = 1.0" << endl;

    matrix2[3][0] = 0.0;
    matrix2[3][1] = 0.0;
    matrix2[3][2] = 378.5;
    matrix2[3][3] = -800.2326;
    matrix2[3][4] = 421.5;
    matrix2[3][5] = 0.0;
    matrix2[3][6] = 0.0;

    cout << "  378.5 * y_2 - 800.2326 * y_3 + 421.5 * y_4 = 1.0" << endl;

    matrix2[4][0] = 0.0;
    matrix2[4][1] = 0.0;
    matrix2[4][2] = 0.0;
    matrix2[4][3] = 378.0;
    matrix2[4][4] = -800.2273;
    matrix2[4][5] = 422.0;
    matrix2[4][6] = 0.0;

    cout << "  378.0 * y_3 - 800.2273 * y_4 + 422.0 * y_5 = 1.0" << endl;

    matrix2[5][0] = 0.0;
    matrix2[5][1] = 0.0;
    matrix2[5][2] = 0.0;
    matrix2[5][3] = 0.0;
    matrix2[5][4] = 377.5;
    matrix2[5][5] = -800.2222;
    matrix2[5][6] = 422.5;

    cout << "  377.5 * y_4 - 800.2222 * y_5 + 422.5 * y_6 = 1.0" << endl;

    matrix2[6][0] = 0.0;
    matrix2[6][1] = 0.0;
    matrix2[6][2] = 0.0;
    matrix2[6][3] = 0.0;
    matrix2[6][4] = 0.0;
    matrix2[6][5] = 0.0;
    matrix2[6][6] = 1.0;

    cout << "  1.0 * y_6 = 2.1599" << endl << endl;

    printMatrix(matrix2, 7);
    cout << endl;

    vector2[0] = -0.05;
    vector2[1] = integral();
    vector2[2] = integral();
    vector2[3] = integral();
    vector2[4] = integral();
    vector2[5] = integral();
    vector2[6] = 2.1599;

    printVector(vector2, 7);
    cout << endl;

    double* y0_6 = solveSystemOfLinearEquations(matrix2, vector2, &cond, 7);

    cout << "  cond = ";
    printf("%10.4f", cond);
    cout << endl << endl;

    cout << "  ";
    for (int i = 0; i < 7; i++) {
        cout << "x_" << i << " = ";
        printf("%4.2f     ", 2.0 + i * 0.05);
    }
    cout << endl;

    cout << "  ";
    for (int i = 0; i < 7; i++) {
        cout << "y_" << i << " = ";
        printf("%6.4f   ", y0_6[i]);
    }
    cout << endl << endl;

    //Решение задачи Коши

    cout << "  RKF45:" << endl << endl;

    int Negn = 2;
    int iwork[30];
    float work[15];
    float Y0[2];
    Y0[0] = y0_3[0];
    Y0[1] = (-y0_3[2] + 4.0 * y0_3[1] - 3.0 * y0_3[0]) / (2.0 * 0.1);
    float T0 = 2.0;
    float RE = 1e-10;
    float AE = 1e-10;
    int iflag = 1;
    float tout = 2.0;
    float h = 0.1;

    float y[4];

    for (int i = 0; i <= 3; i++)
    {
        RKF45(MyFunc, Negn, Y0, &T0, &tout, &RE, &AE, &iflag, work, iwork);
        y[i] = Y0[0];
        tout += h;
    }

    cout << "  ";
    for (int i = 0; i < 4; i++) {
        cout << "x_" << i << " = ";
        printf("%3.1f      ", 2.0 + i * 0.1);
    }
    cout << endl;

    cout << "  ";
    for (int i = 0; i < 4; i++) {
        cout << "y_" << i << " = ";
        printf("%6.4f   ", y[i]);
    }
    cout << endl;
}