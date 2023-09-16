#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define n_anchor 6

typedef struct coordinate{
    double x;
    double y;
    double z;
}coordinate;


double det_3(double a[3][3]){
    double ans = 0;
    ans = a[0][0]*a[1][1]*a[2][2] + a[1][0]*a[2][1]*a[0][2] + a[2][0]*a[0][1]*a[1][2] - a[0][2]*a[1][1]*a[2][0] - a[0][1]*a[1][0]*a[2][2] - a[0][0]*a[1][2]*a[2][1];
    return ans; 
}

double det_2(double *a){
    double ans = 0;
    ans = a[0]*a[3] - a[1]*a[2];
    return ans;
}



double err(coordinate a, coordinate b){
    double ans = sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
    return ans; 
}


int main(){
    time_t t;
    srand((unsigned)time( &t));
    coordinate a0 = {0, 0, 0}, a1 = {5000, 0, 0}, a2 = {0, 5000, 0}, a3 = {5000, 5000, 0}, a4 = {0, 0, 1000}, a5 = {5000, 0, 1000}, a6 = {0, 5000, 1000}, a7 = {5000, 5000, 1000};
    coordinate tag = {1010, 2020, 705};
    coordinate anchor[n_anchor] = {a0, a1, a2, a3, a4, a5};
    double dis[n_anchor] = {0};
    double min = -20;
    double max = 20;
    for(int i = 0; i < n_anchor; i++){
        double error = (max - min) * rand() / (RAND_MAX + 1.0) + min;
        printf("err = %f ", error);
        dis[i] = sqrt(pow(tag.x - anchor[i].x, 2) + pow(tag.y - anchor[i].y, 2) + pow(tag.z - anchor[i].z, 2)) + error;
    }
    printf("\n");
    
    // calculate A
    double A[n_anchor-1][3], b[n_anchor-1], AT[3][n_anchor-1], Ainv[3][3]={0}, Amul[3][3] = {0};
    for(int i = 1; i < n_anchor; i++){
        A[i-1][0] = 2 * (anchor[i].x - anchor[0].x);
        A[i-1][1] = 2 * (anchor[i].y - anchor[0].y);
        A[i-1][2] = 2 * (anchor[i].z - anchor[0].z);
    }
   
    // calculate xi^2 + yi^2 + zi^2
    double k[n_anchor];
    for(int i = 0; i < n_anchor; i++){
        k[i] = pow(anchor[i].x, 2) + pow(anchor[i].y, 2) + pow(anchor[i].z, 2);
    }

    // calculate b 
    for(int i = 0; i < n_anchor-1; i++){
        b[i] = pow(dis[0], 2) - pow(dis[i+1], 2) -  k[0] + k[i+1];
    }
    // calculate AT
   
    for(int i = 0; i < n_anchor - 1; i++){
        for(int j = 0; j < 3; j++){
            AT[j][i] = A[i][j];
        }
    }
    
    //print_3m(AT);
    //calculate Amul = AT * A 
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < n_anchor - 1; k++){
                Amul[i][j] += AT[i][k] * A[k][j];
            }
        }
    }

    //calculate (AT * A)^-1

    double adj[4] = {0};
    double coefficent = 1 / det_3(Amul);
    int index = 0;
    int flag = 0;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                for(int l = 0; l < 3; l++){
                    if(k != i){
                        if(l != j){
                            adj[index] = Amul[k][l];
                            index++;
                        }
                    }
                }
            }
            index = 0;
            if(flag == 0){
                Ainv[i][j] = coefficent * det_2(adj);
                flag = 1;
            }
            else{
                Ainv[i][j] = -(coefficent * det_2(adj));
                flag = 0;
            }
        }
    }   
    
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 3; j++){
            if(i < j){
                double temp = Ainv[i][j];
                Ainv[i][j] = Ainv[j][i];
                Ainv[j][i] = temp;
            }
        }
    }
    //calculate Ainv * AT * b

    double c[3][n_anchor-1] = {0};
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < n_anchor-1; j++){
            for(int k = 0; k < 3; k++){
                c[i][j] += Ainv[i][k]*AT[k][j];
            }
        }
    }
    double cal[3] = {0};
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < n_anchor-1; j++){
            cal[i] += c[i][j] * b[j];
        }
    }
    printf("calculation\n");
    print_mat(cal);
    printf("%f %f %f\n", tag.x, tag.y, tag.z);
    coordinate calculate = {cal[0], cal[1], cal[2]};
    printf("%f %f %f\n", tag.x-calculate.x, tag.y-calculate.y, tag.z-calculate.z);
    //printf("error = %f\n", err(calculate, tag));
    
    return 0;
}

