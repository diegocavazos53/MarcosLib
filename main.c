// Programa para la solución de marcos planos
// incluyendo liberaciones de momentos
// Análisis Matricial de Estructuras
// FIC - UANL - DCL (11-10-2025)

// Directivas del preprocesador
#include <stdio.h>  // Entrada y salida.
#include <stdlib.h> // Librería estándar.
#include <math.h>   // Librería matemática.
#include <string.h> // Manejo de cadenas de texto.

// Constantes
#define ncjt 3 // Número de coordenadas máximas restringidas por nodo.
#define tam 6  // Número de coordenadas por barra.

// Variables globales
double bk[tam][tam];  // Matriz de rigidez de barra en SCL.
double gk[tam][tam];  // Matriz de rigidez de barra en SCG.
double s[1000][1000]; // Matriz de rigidez del sistema.
double qf[tam];       // Vector de cargas de extremo fijo en SCL.
double ff[tam];       // Vector de cargas de extremo fijo en SCG.
double p[1000];       // Vector de cargas del sistema.
double u[tam];        // Vector de desplazamientos en el SCL.
double v[tam];        // Vector de desplazamientos en el SCG.
double ts[tam][tam];  // Matriz de transformación intermedia.
double t[tam][tam];   // Matriz de transformación.
double q[tam];        // Vector de fuerzas en el SCL.
double f[tam];        // Vector de fuerzas en el SCG.
double r[100];        // Vector de reacciones.
double mef[100][6];   // Contenedor de fuerzas en los extremos de los miembros (impr).
int nodoInFin[100][2];// Contenedor de números de nodos (para impresión).

// Función para obtener la matriz de rigidez de barra en el SCG.
void mstiffg(){
    // Inicializar en 0 la matriz de transformación intermedia y la matriz de rigidez
    // de barra en el SCG.
    for(int i = 0; i < 2 * ncjt; i++){
        for(int j = 0; j < 2 * ncjt; j++){
            ts[i][j] = 0;
            gk[i][j] = 0;
        }
    }

    // Matriz de transformación intermedia.
    for(int i = 0; i < 2 * ncjt; i++){
        for(int j = 0; j < 2 * ncjt; j++){
            for(int k = 0; k < 2 * ncjt; k++){
                ts[i][j] = ts[i][j] + bk[i][k] * t[k][j];
            }
        }
    }

    // Matriz de rigidez transformada al Sistema de Coordenadas Globales.
    for(int i = 0; i < 2 * ncjt; i++){
        for(int j = 0; j < 2 * ncjt; j++){
            for(int k = 0; k < 2 * ncjt; k++){
                gk[i][j] = gk[i][j] + t[k][i] * ts[k][j];
            }
        }
    }
}

// Función para almacenar las matrices de rigidez en el SCG en la matriz de
// rigidez ensamblada del sistema.
void stores(int jb, int je, int ndof, int nsc[]){
    int i1 = 0;
    int n1 = 0;
    int n2 = 0;
    for (int i = 0; i < 2 * ncjt; i++){
        if (i < ncjt){
            i1 = (jb - 1) * ncjt + i;
        }
        else{
            i1 = (je - 1) * ncjt + (i - ncjt);
        }
        n1 = nsc[i1];
        if(n1 <= ndof){
            for(int j = 0; j < 2 * ncjt; j++){
                if (j < ncjt){
                    i1 = (jb - 1) * ncjt + j;
                }
                else{
                    i1 = (je - 1) * ncjt + (j - ncjt);
                }
                n2 = nsc[i1];
                if(n2 <= ndof){
                    s[n1 - 1][n2 - 1] = s[n1-1][n2-1]+gk[i][j];
                }
            }
        }
    }
}

// Función para obtener los desplazamientos en el SCG.
void mdispg(int jb, int je, int ndof, int nsc[]){
    int j = 0;
    int n = 0;
    for(int i = 0; i < 2 * ncjt; i++){
        v[i] = 0;
    }
    j = (jb - 1) * ncjt;
    for(int i = 0; i < ncjt; i++){
        j = j + 1;
        n = nsc[j - 1];
        if(n <= ndof){
            v[i] = p[n-1];
        }
    }
    j = (je - 1) * ncjt;
    for(int i = ncjt + 1; i <= 2 * ncjt; i++){
        j = j + 1;
        n = nsc[j - 1];
        if(n <= ndof){
            v[i - 1] = p[n - 1];
        }
    }
}

// Función para fabricar la matriz de transformación.
void mtrans(double cx, double cy){
    for(int i = 0; i < 2 * ncjt; i++){
        for(int j = 0; j < 2 * ncjt; j++){
            t[i][j] = 0;
        }
    }
    t[0][0] = cx; t[1][0] = -cy;
    t[0][1] = cy; t[1][1] =  cx;
    t[3][3] = cx; t[4][3] = -cy;
    t[3][4] = cy; t[4][4] =  cx;
    t[2][2] =  1; t[5][5] =   1;
}

// Función para obtener los desplazamientos en el SCL.
void mdispl(){
    for (int i = 0; i < 2 * ncjt; i++){
        u[i] = 0;
    }
    for (int i = 0; i < 2 * ncjt; i++){
        for (int j = 0; j < 2 * ncjt; j++){
            u[i] = u[i] + t[i][j] * v[j];
        }
    }
}

// Función para obtener la matriz de rigidez de barra en el SCL.
void mstiffl(double modelast, double area, double zi, double bl, int mt){
    for(int i = 0; i < 2 * ncjt; i++){
        for(int j = 0; j < 2 * ncjt; j++){
            bk[i][j] = 0;
        }
    }
    double z = modelast * area / bl;
    bk[0][0] = z; bk[3][0] = -z;
    bk[0][3] = -z; bk[3][3] = z;
    z = modelast * zi / pow(bl,3);

    if (mt == 0){
        bk[1][1] =  12 * z;     bk[2][1] = 6 * bl * z;
        bk[4][1] = -12 * z;     bk[5][1] = 6 * bl * z;
        bk[1][2] =  6 * bl * z; bk[2][2] = 4 * bl * bl * z;
        bk[4][2] = -6 * bl * z; bk[5][2] = 2 * bl * bl * z;
        bk[1][4] = -12 * z;     bk[2][4] = -6 * bl * z;
        bk[4][4] =  12 * z;     bk[5][4] = -6 * bl * z;
        bk[1][5] =  6 * bl * z; bk[2][5] = 2 * bl * bl * z;
        bk[4][5] = -6 * bl * z; bk[5][5] = 4 * bl * bl * z;
    }
    else if (mt == 1){
        bk[1][1] =  3 * z;      bk[4][1] = -3 * z;
        bk[5][1] =  3 * bl * z; bk[1][4] = -3 * z;
        bk[4][4] =  3 * z;      bk[5][4] = -3 * bl * z;
        bk[1][5] =  3 * bl * z; bk[4][5] = -3 * bl * z;
        bk[5][5] =  3 * bl * bl * z;
    }
    else if (mt == 2){
        bk[1][1] =  3 * z;      bk[2][1] =  3 * bl * z;
        bk[4][1] = -3 * z;      bk[1][2] =  3 * bl * z;
        bk[2][2] =  3 * bl * bl * z;
        bk[4][2] = -3 * bl * z; bk[1][4] = -3 * z;
        bk[2][4] = -3 * bl * z; bk[4][4] =  3 * z;
    }
    else if (mt == 3){
        // Ya están los coeficientes de armaduras.
    }
}

// Función para obtener el vector de fuerzas en el SCL.
void mforcel(){
    for (int i = 0; i < 2 * ncjt; i++){
        q[i] = qf[i];
    }
    for (int i = 0; i < 2 * ncjt; i++){
        for (int j = 0; j < 2 * ncjt; j++){
            q[i] = q[i] + bk[i][j] * u[j];
        }
    }
}

// Función para obtener el vector de fuerzas en el SCG.
void mforceg(){
    for (int i = 0; i < 2 * ncjt; i++){
        f[i] = 0;
    }
    for (int i = 0; i < 2 * ncjt; i++){
        for (int j = 0; j < 2 * ncjt; j++){
            f[i] = f[i] + t[j][i] * q[j];
        }
    }
}

// Función para almacenar el vector de reacciones.
void storer(int jb, int je, int ndof, int nsc[]){
    int i1 = 0;
    int n = 0;
    for (int i = 1; i <= 2 * ncjt; i++){
        if(i <= ncjt){
            i1 = (jb - 1) * ncjt + i;
        }
        else{
            i1 = (je - 1) * ncjt + (i - ncjt);
        }
        n = nsc[i1 - 1];
        if(n > ndof){
            r[n - ndof - 1] = r[n - ndof - 1] + f[i - 1];
        }
    }
}

// Función para obtener el vector de fuerzas de extremo fijo en el SCL.
void mfefll(int iml, double bl, int mp[][2], double pm[][4], int mt){
    double fab = 0, fsb = 0, fmb = 0, fae = 0, fse = 0, fme = 0;
    double bw, bl1, bm, w, bl2, w1, w2;
    int ldtype = 0;
    ldtype = mp[iml][1];
    if(ldtype == 1){
        bw = pm[iml][0]; bl1 = pm[iml][2];
        bl2 = bl - bl1;
        fab = 0;
        fsb = bw*pow(bl2,2)/pow(bl,3)*(3*bl1+bl2);
        fmb = bw*bl1*pow(bl2,2)/pow(bl,2);
        fae = 0;
        fse = bw*pow(bl1,2)/pow(bl,3)*(bl1+3*bl2);
        fme = -bw*bl2*pow(bl1,2)/pow(bl,2);
        printf("Fuerza concentrada en el miembro %d:\n", mp[iml][0]);
        printf("P =%13.2f\n", bw);
        printf("L =%13.2f\n", bl);
        printf("l1 =%12.2f\n", bl1);
        printf("l2 =%12.2f\n", bl2);
        printf("fab =%11.2f\n", fab);
        printf("fsb =%11.2f\n", fsb);
        printf("fmb =%11.2f\n", fmb);
        printf("fae =%11.2f\n", fae);
        printf("fse =%11.2f\n", fse);
        printf("fme =%11.2f\n", fme);
        printf("\n");
    }
    if(ldtype == 2){
        bm = pm[iml][0]; bl1 = pm[iml][2];
        bl2 = bl - bl1;
        fab = 0;
        fsb = -6*bm*bl1*bl2/pow(bl,3);
        fmb = bm*bl2/pow(bl,2)*(bl2-2*bl1);
        fae = 0;
        fse = 6*bm*bl1*bl2/pow(bl,3);
        fme = bm*bl1/pow(bl,2)*(bl1-2*bl2);
        printf("Momento concentrado en el miembro %d:\n", mp[iml][0]);
        printf("M =%13.2f\n", bm);
        printf("L =%13.2f\n", bl);
        printf("l1 =%12.2f\n", bl1);
        printf("l2 =%12.2f\n", bl2);
        printf("fab =%11.2f\n", fab);
        printf("fsb =%11.2f\n", fsb);
        printf("fmb =%11.2f\n", fmb);
        printf("fae =%11.2f\n", fae);
        printf("fse =%11.2f\n", fse);
        printf("fme =%11.2f\n", fme);
        printf("\n");
    }
    if(ldtype == 3){
        w = pm[iml][0]; bl1 = pm[iml][2]; bl2 = pm[iml][3];
        fab = 0;
        fsb = w*bl/2*(1-bl1/pow(bl,4)*(2*pow(bl,3)-2*pow(bl1,2)*bl+pow(bl1,3))-pow(bl2,3)/pow(bl,4)*(2*bl-bl2));
        fmb = w*pow(bl,2)/12*(1-pow(bl1,2)/pow(bl,4)*(6*pow(bl,2)-8*bl1*bl+3*pow(bl1,2))-pow(bl2,3)/pow(bl,4)*(4*bl-3*bl2));
        fae = 0;
        fse = w*bl/2*(1-pow(bl1,3)/pow(bl,4)*(2*bl-bl1)-bl2/pow(bl,4)*(2*pow(bl,3)-2*pow(bl2,2)*bl+pow(bl2,3)));
        fme = -w*pow(bl,2)/12*(1-pow(bl1,3)/pow(bl,4)*(4*bl-3*bl1)-pow(bl2,2)/pow(bl,4)*(6*pow(bl,2)-8*bl2*bl+3*pow(bl2,2)));
        printf("Carga uniformemente repartida en el miembro %d:\n", mp[iml][0]);
        printf("w =%13.2f\n", w);
        printf("L =%13.2f\n", bl);
        printf("l1 =%12.2f\n", bl1);
        printf("l2 =%12.2f\n", bl2);
        printf("fab =%11.2f\n", fab);
        printf("fsb =%11.2f\n", fsb);
        printf("fmb =%11.2f\n", fmb);
        printf("fae =%11.2f\n", fae);
        printf("fse =%11.2f\n", fse);
        printf("fme =%11.2f\n", fme);
        printf("\n");
    }
    if(ldtype == 4){
        w1 = pm[iml][0]; w2 = pm[iml][1]; bl1 = pm[iml][2]; bl2 = pm[iml][3];
        fab = 0;
        fsb = w1*pow(bl-bl1,3)/20/pow(bl,3)*((7*bl+8*bl1)-bl2*(3*bl+2*bl1)/(bl-bl1)*(1+bl2/(bl-bl1)+pow(bl2,2)/pow(bl-bl1,2))+2*pow(bl2,4)/pow(bl-bl1,3));
        fsb = fsb + w2*pow(bl-bl1,3)/20/pow(bl,3)*((3*bl+2*bl1)*(1+bl2/(bl-bl1)+pow(bl2,2)/pow(bl-bl1,2))-pow(bl2,3)/pow(bl-bl1,2)*(2+(15*bl-8*bl2)/(bl-bl1)));
        fmb = w1*pow(bl-bl1,3)/60/pow(bl,2)*((3*bl+4*bl1)-bl2*(2*bl+3*bl1)/(bl-bl1)*(1+bl2/(bl-bl1)+pow(bl2,2)/pow(bl-bl1,2))+3*pow(bl2,4)/pow(bl-bl1,3));
        fmb = fmb + w2*pow(bl-bl1,3)/60/pow(bl,2)*((2*bl+3*bl1)*(1+bl2/(bl-bl1)+pow(bl2,2)/pow(bl-bl1,2))-3*pow(bl2,3)/pow(bl-bl1,2)*(1+(5*bl-4*bl2)/(bl-bl1)));
        fae = 0;
        fse = (w1+w2)/2*(bl-bl1-bl2)-fsb;
        fme = (bl-bl1-bl2)/6*(w1*(-2*bl+2*bl1-bl2)-w2*(bl-bl1+2*bl2))+fsb*bl-fmb;
        printf("Carga variablemente repartida en el miembro %d:\n", mp[iml][0]);
        printf("w1 =%12.2f\n", w1);
        printf("w2 =%12.2f\n", w2);
        printf("L =%13.2f\n", bl);
        printf("l1 =%12.2f\n", bl1);
        printf("l2 =%12.2f\n", bl2);
        printf("fab =%11.2f\n", fab);
        printf("fsb =%11.2f\n", fsb);
        printf("fmb =%11.2f\n", fmb);
        printf("fae =%11.2f\n", fae);
        printf("fse =%11.2f\n", fse);
        printf("fme =%11.2f\n", fme);
        printf("\n");
    }
    if(ldtype == 5){
        bw = pm[iml][0]; bl1 = pm[iml][2];
        bl2 = bl - bl1;
        fab = bw*bl2/bl;
        fsb = 0;
        fmb = 0;
        fae = bw*bl1/bl;
        fse = 0;
        fme = 0;
        printf("Fuerza axial concentrada en el miembro %d:\n", mp[iml][0]);
        printf("P =%13.2f\n", bw);
        printf("L =%13.2f\n", bl);
        printf("l1 =%12.2f\n", bl1);
        printf("l2 =%12.2f\n", bl2);
        printf("fab =%11.2f\n", fab);
        printf("fsb =%11.2f\n", fsb);
        printf("fmb =%11.2f\n", fmb);
        printf("fae =%11.2f\n", fae);
        printf("fse =%11.2f\n", fse);
        printf("fme =%11.2f\n", fme);
        printf("\n");
    }
    if(ldtype == 6){
        w = pm[iml][0]; bl1 = pm[iml][2]; bl2 = pm[iml][3];
        fab = w/2/bl*(bl-bl1-bl2)*(bl-bl1+bl2);
        fsb = 0;
        fmb = 0;
        fae = w/2/bl*(bl-bl1-bl2)*(bl+bl1-bl2);
        fse = 0;
        fme = 0;
        printf("Fuerza axial repartida en el miembro %d:\n", mp[iml][0]);
        printf("w =%13.2f\n", w);
        printf("L =%13.2f\n", bl);
        printf("l1 =%12.2f\n", bl1);
        printf("l2 =%12.2f\n", bl2);
        printf("fab =%11.2f\n", fab);
        printf("fsb =%11.2f\n", fsb);
        printf("fmb =%11.2f\n", fmb);
        printf("fae =%11.2f\n", fae);
        printf("fse =%11.2f\n", fse);
        printf("fme =%11.2f\n", fme);
        printf("\n");
    }
    qf[0] = qf[0] + fab;
    qf[3] = qf[3] + fae;

    if (mt == 0){
        qf[1] = qf[1] + fsb;
        qf[2] = qf[2] + fmb;
        qf[4] = qf[4] + fse;
        qf[5] = qf[5] + fme;
    }
    else if (mt == 1){
        qf[1] = qf[1] + fsb - 3.0 * fmb/(2.0*bl);
        qf[4] = qf[4] + fse + 3.0 * fmb/(2.0*bl);
        qf[5] = qf[5] + fme - fmb/2.0;
    }
    else if (mt == 2){
        qf[1] = qf[1] + fsb - 3.0 * fme/(2.0*bl);
        qf[2] = qf[2] + fmb - fme/2.0;
        qf[4] = qf[4] + fse + 3.0 * fme/(2.0*bl);
    }
    else if (mt == 3){
        qf[1] = qf[1] + fsb - (fmb+fme)/bl;
        qf[4] = qf[4] + fse + (fmb+fme)/bl;
    }
}

// Función para obtener el vector de fuerzas de extremo fijo en el SCG.
void mfefg(){
    for(int i = 0; i < 2 * ncjt; i++){
        ff[i] = 0;
    }
    for(int i = 0; i < 2 * ncjt; i++){
        for(int j = 0; j < 2 * ncjt; j++){
            ff[i] = ff[i] + t[j][i] * qf[j];
        }
    }
}

// Función para ensamblar el vector de fuerzas de extremo fijo del sistema.
void storepf(int jb, int je, int ndof, int nsc[]){
    int i1 = 0;
    int n1 = 0;
    for(int i = 0; i < 2 * ncjt; i++){
        if(i < ncjt){
            i1 = (jb - 1) * ncjt + i;
        }
        else{
            i1 = (je - 1) * ncjt + (i - ncjt);
        }
        n1 = nsc[i1];
        if(n1 <= ndof){
            p[n1 - 1] = p[n1 - 1] - ff[i];
        }
    }
}

// Programa principal
/* Se le pueden pasar argumentos desde la línea de comandos,
   como por ejemplo, el nombre de archivo que se desea leer.
*/
int main(int argc, char** argv)
{
    printf("===  M A R C O S   P L A N O S  ===\n");
    printf("incluyendo liberaciones de momentos\n\n");
    printf("===  Diego Cavazos de Lira 2025 ===\n\n");
    printf("Argumentos: %d\n",argc); // Número de argumentos en la consola.
    // El argumento [0] es la dirección del ejecutable.
    // El argumento [1] lo utilizaremos para indicar un archivo a leer.
    printf("Argumentos: %s\n\n",argv[1]);
    char nomArchivo[50];
    // Si no se especifica desde la terminal qué archivo leer, preguntar
    // al usuario el nombre del archivo a leer.
    if (argc > 1){
        strcpy(nomArchivo, argv[1]);
    }
    else{
        printf("¿Qué archivo desea leer?\n");
        scanf("%s", nomArchivo);
    }

    // Si se desea siempre leer el mismo archivo, su nombre se puede poner
    // a continuación.
    // char nomArchivo[25] = "ma1.m2d";
    printf("Leyendo: %s\n\n",nomArchivo);

    FILE *pFile;
    pFile = fopen(nomArchivo,"r");
    if(!pFile){
        printf("ERROR: No se pudo leer el archivo.\n");
        return 1;
    }

    // Parte 1
    int nj = 0; // Número de nodos.
    fscanf(pFile,"%d", &nj);

    double coord[nj][2]; // Matriz de coordenadas.
    for (int i = 0; i < nj; i++){
        fscanf(pFile,"%lf",&coord[i][0]);
        fscanf(pFile,"%lf",&coord[i][1]);
    }

    // Parte 2
    int ns = 0; // Número de apoyos.
    fscanf(pFile,"%d", &ns);

    int msup[ns][ncjt + 1]; // Matriz de apoyos.
    for (int i = 0; i < ns; i++){
        for(int j = 0; j < ncjt + 1; j++){
            fscanf(pFile,"%d",&msup[i][j]);
        }
    }

    // Parte 3
    int nmp = 0; // Número de propiedades de los materiales.
    fscanf(pFile,"%d", &nmp);

    double em[nmp]; // Vector de módulos elásticos.
    for (int i = 0; i < nmp; i++){
        fscanf(pFile,"%lf",&em[i]);
    }

    // Parte 4
    int ncp = 0; // Número de secciones transversales.
    fscanf(pFile,"%d", &ncp);

    double cp[ncp][2]; // Ahora se tiene el área y la inercia.
    for (int i = 0; i < ncp; i++){
        fscanf(pFile,"%lf",&cp[i][0]); // Área.
        fscanf(pFile,"%lf",&cp[i][1]); // Inercia.
    }

    // Parte 5
    int nm = 0; // Número de miembros.
    fscanf(pFile,"%d", &nm);

    int mprp[nm][5]; // Matriz de propiedades de los miembros.
    for (int i = 0; i < nm; i++){
        for(int j = 0; j < 5; j++){
            fscanf(pFile,"%d",&mprp[i][j]);
        }
    }

    // Parte 6
    int njl = 0; // Número de cargas nodales.
    fscanf(pFile,"%d", &njl);

    int jp[njl]; // Vector con los números de los nodos cargados.
    double pj[njl][ncjt]; // Cargas en los nodos (Fx, Fy, Mz).
    for(int i = 0; i < njl; i++){
        fscanf(pFile,"%d", &jp[i]);
        for(int j = 0; j < ncjt; j++){
            fscanf(pFile,"%lf",&pj[i][j]);
        }
    }

    // Parte 6B
    int nml = 0; // Número de cargas en de los miembros.
    fscanf(pFile,"%d", &nml);

    int mp[nml][2]; // Vector con los números de los miembros cargados y tipo de carga.
    double pm[nml][4]; // Cargas en los miembros.
    for(int i = 0; i < nml; i++){
        for(int j = 0; j < 2; j++){
            fscanf(pFile,"%d",&mp[i][j]);
        }
        for(int j = 2; j < 6; j++){
            fscanf(pFile,"%lf",&pm[i][j-2]);
        }
    }

    printf("El archivo se leyó correctamente.\n");

    if(fclose(pFile) != 0){
        printf("Error: No se pudo cerrar el archivo.\n");
    }
    printf("\n");

    // Impresión de la información de entrada
    printf("     **********************************\n");
    printf("     *    Programa de Marcos Planos   *\n");
    printf("     *              por               *\n");
    printf("     *      Diego Cavazos de Lira     *\n");
    printf("     **********************************\n\n");
    printf("        =============================\n");
    printf("        DATOS ESTRUCTURALES GENERALES\n");
    printf("        =============================\n\n");

    printf("Tipo de estructura: Marco Plano\n");
    printf("Número de nodos: %d\n", nj);
    printf("Número de miembros: %d\n", nm);
    printf("Número de conjuntos de propiedades de materiales: %d\n", nmp);
    printf("Número de conjuntos de propiedades de sección transversal: %d\n\n", ncp);

    printf("        =============================\n");
    printf("          Coordenadas de los nodos\n");
    printf("        =============================\n");
    printf("No. Nodo     Coordenada X         Coordenada Y\n");
    printf("--------     ------------         ------------\n");
    for(int i = 0; i < nj; i++){
        printf("%4d %20.4lf %20.4lf\n", i+1, coord[i][0], coord[i][1]);
    }
    printf("\n");
    printf("        =============================\n");
    printf("                   Apoyos\n");
    printf("        =============================\n");
    printf("No. Nodo      Restricción X    Restricción Y    Restr. Rotacional\n");
    printf("--------      -------------    -------------    -----------------\n");
    for(int i = 0; i < ns; i++){
        printf("%4d %16d %16d %18d\n", msup[i][0], msup[i][1], msup[i][2], msup[i][3]);
    }
    printf("\n");
    printf("        =============================\n");
    printf("        Propiedades de los materiales\n");
    printf("        =============================\n");
    printf("Número        Módulo de       Coeficiente de\n"); // El coef. térmico se implementará luego.
    printf("material     Elasticidad     expansión térmica\n");
    printf("--------     -----------     -----------------\n");
    for(int i = 0; i < nmp; i++){
        printf("%4d %19.4lf\n", i+1, em[i]);
    }
    printf("\n");
    printf("        =============================\n");
    printf("        Propiedades de secc. transv.\n");
    printf("        =============================\n");
    printf("No. Prop.      Area (A)     Momento de inercia (I)\n");
    printf("---------     ---------     ----------------------\n");
    for(int i = 0; i < ncp; i++){
        printf("%4d %18.4lf %26.4lf\n", i+1, cp[i][0], cp[i][1]);
    }
    printf("\n");
    printf("        =============================\n");
    printf("         Información de los miembros\n");
    printf("        =============================\n");
    printf("Miembro  Nodo de  Nodo de  Número    Número Prop.   Tipo de\n");
    printf("número   Inicio      Fin   material  secc. transv.  miembro\n");
    printf("-------  -------  -------  --------  -------------  -------\n");
    for(int i = 0; i < nm; i++){
        printf("%4d %8d %8d %9d %11d %11d\n", i+1, mprp[i][0], mprp[i][1], mprp[i][2], mprp[i][3], mprp[i][4]);
    }
    printf("\n");
    printf("        =============================\n");
    printf("                Cargas nodales\n");
    printf("        =============================\n");
    printf("No. Nodo       Fuerza X           Fuerza Y            Momento\n");
    printf("--------       --------           --------            -------\n");
    for(int i = 0; i <njl; i++){
        printf("%4d %18.4lf %18.4lf %18.2lf\n", jp[i], pj[i][0], pj[i][1], pj[i][2]);
    }
    printf("\n");
    printf("        =============================\n");
    printf("            Cargas en los miembros\n");
    printf("        =============================\n");
    printf("                                        Magnitud\n");
    printf("                                      de la carga\n");
    printf("                                       (W o M) o     Intensidad\n");
    printf("                     No.     ID de     Intensidad   de la carga   Distancia   Distancia\n");
    printf(" Tipo de carga    Miembro    carga      (w o w1)         w2           l1          l2\n");
    printf("---------------   -------   -------   -----------   -----------   ---------   ---------\n");
    for(int i = 0; i <nml; i++){
        if(mp[i][1] == 1){
            printf("F. Puntual     ");
        }
        else if(mp[i][1] == 2){
            printf("M. Puntual     ");
        }
        else if(mp[i][1] == 3){
            printf("Unif. Repartida");
        }
        else if(mp[i][1] == 4){
            printf("Var. Repartida ");
        }
        else if(mp[i][1] == 5){
            printf("F. Axial       ");
        }
        else if(mp[i][1] == 6){
            printf("F. Axial Rep.  ");
        }
        else if(mp[i][1] == 7){
            printf("M. Torsionante ");
        }
        else{
            printf("¿?             ");
        }
        printf("%7d %9d %16.4lf %13.4lf %11.4lf %11.4lf\n", mp[i][0], mp[i][1], pm[i][0], pm[i][1], pm[i][2], pm[i][3]); // Revisar centrado.
    }
    printf("\n");
    printf("****** Fin de la información de entrada ******\n\n");

    // Módulo de Análisis Estructural
    // Parte 7
    int nr = 0;
    int ndof = 0;
    for (int i = 0; i < ns; i++){
        for (int j = 1; j < ncjt + 1; j++){
            if (msup[i][j] == 1){
                nr = nr + 1;
            }
        }
    }
    ndof = ncjt * nj - nr;
    printf("El número de coordenadas restringidas NR es: %d\n", nr);
    printf("El número de grados de libertad NDOF es: %d\n", ndof);

    // Parte 8
    int nsc[ncjt*nj];
    int icount = 0;
    int i3 = 0;
    int j = 0;
    int k = ndof;
    for(int i = 1; i <= nj; i++){
        icount = 0;
        for(int i1 = 0; i1 < ns; i1++){
            if(msup[i1][0] == i){
                icount = 1;
                for(int i2 = 0; i2 < ncjt; i2++){
                    i3 = (i - 1) * ncjt + i2;
                    if(msup[i1][i2+1] == 1){
                        k = k + 1;
                        nsc[i3] = k;
                    }
                    else{
                        j = j + 1;
                        nsc[i3] = j;
                    }
                }
            }
        }
        if (icount == 0){
            for(int i2 = 0; i2 < ncjt; i2++){
                i3 = (i-1) * ncjt + i2;
                j = j + 1;
                nsc[i3] = j;
            }
        }
    }
    printf("El vector de números de coordenadas de la estructura NSC, es:\n");
    for(int i = 0; i < ncjt * nj; i++){
        printf("%4d\n", nsc[i]);
    }
    printf("\n");

    // Parte 9
    int jb = 0; int je = 0;
    int imatnum = 0; int isecctransvnum = 0;
    double modelast = 0;
    double area = 0; double zi = 0; // Inercia Izz.
    double xb = 0; double xe = 0;
    double yb = 0; double ye = 0;
    double bl = 0;
    double cx = 0; double cy = 0;
    int mt = 0;

    for (int im = 0; im < nm; im++){
        jb = mprp[im][0];
        je = mprp[im][1];
        imatnum = mprp[im][2];
        isecctransvnum = mprp[im][3];
        mt = mprp[im][4];
        modelast = em[imatnum - 1];
        area = cp[isecctransvnum - 1][0];
        zi = cp[isecctransvnum - 1][1];
        xb = coord[jb - 1][0];
        yb = coord[jb - 1][1];
        xe = coord[je - 1][0];
        ye = coord[je - 1][1];
        bl = sqrt(pow(xe-xb,2)+pow(ye-yb,2));
        cx = (xe-xb)/bl;
        cy = (ye-yb)/bl;

        // Matriz de rigidez en el SCL.
        mstiffl(modelast, area, zi, bl, mt);
        printf("La matriz de rigidez bk de la barra no. %d es:\n",im + 1);
        for(int i = 0; i < 2 * ncjt; i++){
            for(int j = 0; j < 2 * ncjt; j++){
                printf("%15.3f ", bk[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        // Matriz de transformación.
        mtrans(cx, cy);
        printf("La matriz de transformación t de la barra no. %d es:\n",im + 1);
        for(int i = 0; i < 2 * ncjt; i++){
            for(int j = 0; j < 2 * ncjt; j++){
                printf("%15.3f ", t[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        // Matriz de rigidez en el SCG.
        mstiffg();
        printf("La matriz de rigidez en el SCG de la barra %d es:\n", im + 1);
        for(int i = 0; i < 2*ncjt; i++){
            for(int j = 0; j < 2* ncjt; j++){
                printf("%15.3f ", gk[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        // Almacenar la matriz de rigidez del sistema.
        stores(jb, je, ndof, nsc);

        // Atención a cargas en los miembros.
        if(nml > 0){
            for(int i = 0; i < 2 * ncjt; i++){
                qf[i] = 0;
            }
            for(int iml = 0; iml < nml; iml++){
                if(im + 1 == mp[iml][0]){
                    mfefll(iml, bl, mp, pm, mt);
                }
            }
            printf("El vector de fuerzas locales Qf de la barra no. %d es:\n", im + 1);
            for(int i = 0; i < 2 * ncjt; i++){
                printf("%12.2f\n", qf[i]);
            }
            printf("\n");

            mfefg();
            printf("El vector de fuerzas globales Ff de la barra no. %d es:\n", im + 1);
            for(int i = 0; i < 2 * ncjt; i++){
                printf("%12.2f\n", ff[i]);
            }
            printf("\n");

            storepf(jb, je, ndof, nsc);
        }
    }
    printf("La matriz de rigidez del sistema es:\n");
    for(int i = 0; i < ndof; i++){
        for(int j = 0; j < ndof; j++){
            printf("%15.3f ", s[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("El negativo del vector de fuerzas de extremo fijo en el SCG, -Pf, es:\n");
    for(int i = 0; i < ndof; i++){
        printf("%15.3f\n", p[i]);
    }
    printf("\n");

    // Parte 10
    for (int i = 0; i < njl; i++){
        int i1 = jp[i];
        int i2 = (i1-1) * ncjt;
        int n = 0;
        for (int j = 0; j < ncjt; j++){
            i2 = i2 + 1;
            n = nsc[i2 - 1];
            if (n <= ndof){
                p[n - 1] = p[n - 1] + pj[i][j];
            }
        }
    }

    printf("El vector de fuerzas P-Pf, es:\n");
    for (int i = 0; i < ndof; i++){
        printf("%15.3lf\n", p[i]);
    }
    printf("\n");

    // Parte 11
    double z1 = 0;
    double z = 0;
    for(int i = 0; i < ndof; i++){
        z1 = s[i][i];
        for(int j = 0; j < ndof; j++){
            s[i][j] = s[i][j] / z1;
        }
        p[i] = p[i] / z1;
        for(int k = 0; k < ndof; k++){
            if(k == i){
            }
            else{
                z = s[k][i];
                for(int m = i; m < ndof; m++){
                    s[k][m] = s[k][m] - s[i][m] * z;
                }
                p[k] = p[k] - p[i] * z;
            }
        }
    }

    printf("El vector de desplazamientos d, es:\n");
    for (int i = 0; i < ndof; i++){
        printf("%15.3lf\n", p[i]);
    }
    printf("\n");

    // Parte 12
    for (int im = 0; im < nm; im++){
        jb = mprp[im][0];
        je = mprp[im][1];
        imatnum = mprp[im][2];
        isecctransvnum = mprp[im][3];
        mt = mprp[im][4];
        modelast = em[imatnum - 1];
        area = cp[isecctransvnum - 1][0];
        zi = cp[isecctransvnum - 1][1];
        xb = coord[jb - 1][0];
        yb = coord[jb - 1][1];
        xe = coord[je - 1][0];
        ye = coord[je - 1][1];
        bl = sqrt(pow(xe-xb,2)+pow(ye-yb,2));
        cx = (xe-xb)/bl;
        cy = (ye-yb)/bl;

        // Vector de desplazamientos en el SCG de cada barra.
        mdispg(jb, je, ndof, nsc);
        printf("El vector de desplazamientos globales V de la barra no. %d, es:\n",im + 1);
        for (int i = 0; i < 2 * ncjt; i++){
            printf("%15.3lf\n", v[i]);
        }
        printf("\n");

        mtrans(cx, cy);
        printf("La matriz de transformación de la barra %d es:\n",im+1);
        for (int i = 0; i < 2 * ncjt; i++){
            for (int j = 0; j < 2 * ncjt; j++){
                printf("%15.3lf", t[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        // Vector de desplazamientos en el SCL de cada barra.
        mdispl();
        printf("El vector de desplazamientos en SCL, u de la barra %d es:\n", im+1);
        for (int i = 0; i < 2 * ncjt; i++){
            printf("%15.3lf\n",u[i]);
        }
        printf("\n");

        mstiffl(modelast, area, zi, bl, mt);

        for(int i = 0; i < 2 * ncjt; i++){
            qf[i] = 0;
        }
        if(nml > 0){
            for(int iml = 0; iml < nml; iml++){
                if(im + 1 == mp[iml][0]){
                    mfefll(iml, bl, mp, pm, mt);
                }
            }
        }

        mforcel();
        printf("El vector de fuerzas en el SCL Q de la barra no. %d es:\n", im+1);
        for(int i = 0; i < 2 * ncjt; i++){
            printf("%15.3lf\n", q[i]);
        }
        printf("\n");

        mforceg();
        printf("El vector de fuerzas en el SCG F de la barra no. %d es:\n", im+1);
        for(int i = 0; i < 2 * ncjt; i++){
            printf("%15.3lf\n", f[i]);
        }
        printf("\n");

        storer(jb, je, ndof, nsc);
        for (int i = 0; i < 2 * ncjt; i++){
            mef[im][i] = q[i];
        }
        for (int i = 0; i < 2; i++){
            nodoInFin[im][0] = jb;
            nodoInFin[im][1] = je;
        }
    }

    printf("El vector de reacciones en los apoyos R es:\n");
    for(int i = 0; i < nr; i++){
        printf("%15.3lf\n", r[i]);
    }
    printf("\n");

    // Impresión de salida
    printf("***************************************************\n");
    printf("*             Resultados del Análisis             *\n");
    printf("***************************************************\n\n");

    printf("============================\n");
    printf("Desplazamientos de los nodos\n");
    printf("============================\n\n");

    printf("No. Nodo     Traslación X    Traslación Y    Rotación (rad)\n");
    printf("========     ============    ============    ==============\n");
    double vecdespl[nj*ncjt];
    int contador = 0;
    for(int i = 0; i < nj * ncjt + 1; i++){
        if(nsc[i] > ndof){
            vecdespl[i] = 0;
        }
        if(nsc[i] <= ndof){
            vecdespl[i] = p[contador];
            contador = contador + 1;
        }
    }

    int nnodo = 1;
    for(int i = 0; i < nj * ncjt; i++){
        printf("%4d %20.4lf %15.4lf %17.4lf\n",nnodo,vecdespl[i],vecdespl[i+1], vecdespl[i+2]);
        i = i + 2;
        nnodo = nnodo + 1;
    }
    printf("\n");

    printf("=================================================\n");
    printf("Fuerzas en los Extremos de los Miembros en el SCL\n");
    printf("=================================================\n\n");

    printf("Miembro   Nodo    Fuerza Axial   Fuerza Cortante        Momento\n");
    printf("=======   ====    ============   ===============        =======\n");

    for (int i = 0; i < nm; i++){
        printf("%4d %7d %17.3lf %17.3lf %14.3lf\n", i+1, nodoInFin[i][0], mef[i][0], mef[i][1], mef[i][2]);
        printf("%4d %7d %17.3lf %17.3lf %14.3lf\n\n", i+1, nodoInFin[i][1], mef[i][3], mef[i][4], mef[i][5]);
    }
    printf("\n");

    printf("==============================\n");
    printf("   Reacciones en los apoyos   \n");
    printf("==============================\n\n");

    printf("No. Nodo     Fuerza X       Fuerza Y         Momento\n");
    printf("========     ========       ========         =======\n");
    double vecreacc[ns*ncjt];
    int contadorreacc = 0;
    int contadorceros = 0;
    for (int i = 0; i < ns; i++){
        if(msup[i][1] == 1){
            vecreacc[contadorreacc + contadorceros] = r[contadorreacc];
            contadorreacc = contadorreacc + 1;
        }
        else{
            vecreacc[contadorreacc + contadorceros] = 0;
            contadorceros = contadorceros + 1;
        }
        if(msup[i][2] == 1){
            vecreacc[contadorreacc + contadorceros] = r[contadorreacc];
            contadorreacc = contadorreacc + 1;
        }
        else{
            vecreacc[contadorreacc + contadorceros] = 0;
            contadorceros = contadorceros + 1;
        }
        if(msup[i][3] == 1){
            vecreacc[contadorreacc + contadorceros] = r[contadorreacc];
            contadorreacc = contadorreacc + 1;
        }
        else{
            vecreacc[contadorreacc + contadorceros] = 0;
            contadorceros = contadorceros + 1;
        }
    }

    int indnodo = 0;
    for(int i = 0; i < ns * ncjt; i++){
        printf("%5d %15.3lf %14.3lf %15.3lf\n",msup[indnodo][0],vecreacc[i],vecreacc[i+1],vecreacc[i+2]);
        i = i + 2;
        indnodo = indnodo + 1;
    }
    printf("\n");

    printf("********** Fin del Análisis **********\n");
    getchar();

    return 0;
}
