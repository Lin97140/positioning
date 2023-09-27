#include <math.h>
#include "position.h"

#define N_ANCHOR 4

coordinate anchor[N_ANCHOR] = {{0, 0, 0}, {10, 0, 0}, {0, 10, 0}, {0, 0, 3}};

double determinant_3(double matrix_3[][3]){
    
    return matrix_3[0][0]*matrix_3[1][1]*matrix_3[2][2] + matrix_3[1][0]*matrix_3[2][1]*matrix_3[0][2] + matrix_3[2][0]*matrix_3[0][1]*matrix_3[1][2] - matrix_3[0][2]*matrix_3[1][1]*matrix_3[2][0] - matrix_3[0][1]*matrix_3[1][0]*matrix_3[2][2] - matrix_3[0][0]*matrix_3[1][2]*matrix_3[2][1];
    
}
double determinant_2(double* matrix_2){
    
    return  matrix_2[0]*matrix_2[3] - matrix_2[1]*matrix_2[2];

}

coordinate localization(double distance_matrix[]){

    // calculate A
    double A[N_ANCHOR-1][3], b[N_ANCHOR-1], AT[3][N_ANCHOR-1], ATAinv[3][3]={0}, ATmulA[3][3] = {0};
    for(int i = 1; i < N_ANCHOR; i++){
        A[i-1][0] = 2 * (anchor[i].x - anchor[0].x);
        A[i-1][1] = 2 * (anchor[i].y - anchor[0].y);
        A[i-1][2] = 2 * (anchor[i].z - anchor[0].z);
    }
   
    // calculate xi^2 + yi^2 + zi^2
    double k[N_ANCHOR];
    for(int i = 0; i < N_ANCHOR; i++){
        k[i] = pow(anchor[i].x, 2) + pow(anchor[i].y, 2) + pow(anchor[i].z, 2);
    }

    // calculate b 
    for(int i = 0; i < N_ANCHOR-1; i++){
        b[i] = pow(distance_matrix[0], 2) - pow(distance_matrix[i+1], 2) -  k[0] + k[i+1];
    }

    // calculate AT
    for(int i = 0; i < N_ANCHOR - 1; i++){
        for(int j = 0; j < 3; j++){
            AT[j][i] = A[i][j];
        }
    }
    //calculate calculation = (AT*A)-1*AT*b

    //calculate ATmulA = AT * A 
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < N_ANCHOR - 1; k++){
                ATmulA[i][j] += AT[i][k] * A[k][j];
            }
        }
    }

    //calculate ATAinv = (AT * A)^-1
    double adj[4] = {0};
    double coefficent = 1 / determinant_3(ATmulA);
    int index = 0;
    int flag = 0;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                for(int l = 0; l < 3; l++){
                    if(k != i){
                        if(l != j){
                            adj[index] = ATmulA[k][l];
                            index++;
                        }
                    }
                }
            }
            index = 0;
            if(flag == 0){
                ATAinv[i][j] = coefficent * determinant_2(adj);
                flag = 1;
            }
            else{
                ATAinv[i][j] = -(coefficent * determinant_2(adj));
                flag = 0;
            }
        }
    }   
    
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 3; j++){
            if(i < j){
                double temp = ATAinv[i][j];
                ATAinv[i][j] = ATAinv[j][i];
                ATAinv[j][i] = temp;
            }
        }
    }

    //calculate ATAinv * AT * b
    double c[3][N_ANCHOR-1] = {0};
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < N_ANCHOR-1; j++){
            for(int k = 0; k < 3; k++){
                c[i][j] += ATAinv[i][k]*AT[k][j];
            }
        }
    }
    double calculation[3] = {0};
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < N_ANCHOR-1; j++){
            calculation[i] += c[i][j]*b[j];
        }
    }
    
    //get position
    coordinate position = {calculation[0], calculation[1], calculation[2]};
    return position;
}

