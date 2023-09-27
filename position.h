
typedef struct coordinate{
    double x;
    double y;
    double z;
}coordinate;

double determinant_3(double matrix[][3]);
double determinant_2(double*);

coordinate localization(double distance_matrix[]);
