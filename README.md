# Practica-7-TFIII

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 24 //lado del sistema
#define lado_max 10 //lado máximo de los loops de wilson
#define NormRANu (2.3283063671E-10F)

///Predeclaración de funciones

//Numeros random
void ini_ran(int SEMILLA);
double Random(void);

//Configuracion inicial
void direccionarios(int *xp, int *yp, int *zp, int *xn, int *yn, int *zn);
void configuracion_inicial (int flag, int **spines);

//Configuracion en el equilibrio
int energia_plaqueta_plano_xy (int *xp, int *yp, int **spines, int n, int i, int j);
int energia_plaqueta_plano_xz (int *xp, int *zp, int **spines, int n, int i, int k);
int energia_plaqueta_plano_yz (int *yp, int *zp, int **spines, int n, int j, int k);
int energia_total (int *xp, int *yp, int *zp, int **spines);
void gauge (int *xn, int *yn, int *zn, int **spines, int n);
void Metropolis (int *xp, int *yp, int *zp, int *xn, int *yn, int *zn, int **spines, double *probabilidades, int *energia);

//Wilson
int loop_wilson_xy(int *xp, int *yp, int *xn, int *yn, int **spines, int n, int i, int j, int lado);
int loop_wilson_xz(int *xp, int *zp, int *xn, int *zn, int **spines, int n, int i, int k, int lado);
int loop_wilson_yz(int *xp, int *yp, int *zp, int *xn, int *yn, int *zn, int **spines, int n, int j, int k, int lado);
void estadistica_wilson(int *xp, int *yp, int *zp, int *xn, int *yn, int *zn, int **spines, double *promedio_wilson);


///Números aleatorios, los pongo antes del main porque dependen de variables globales
//Esta funciones las guardo de cuando fisica computacional, ahora solo Dios sabe como chuchas funciona esta cosa
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran, ig1, ig2, ig3;

void ini_ran(int SEMILLA){
    int INI, FACTOR, SUM, i;
    srand(SEMILLA);
    INI = SEMILLA;
    FACTOR = 67397;
    SUM = 7364893;
    for (i = 0; i < 256; i++)
    {
        INI = (INI * FACTOR + SUM);
        irr[i] = INI;
    }
    ind_ran = ig1 = ig2 = ig3 = 0;
}

double Random(void){
    double r;
    ig1 = ind_ran - 24;
    ig2 = ind_ran - 55;
    ig3 = ind_ran - 61;
    irr[ind_ran] = irr[ig1] + irr[ig2];
    ir1 = (irr[ind_ran] ^ irr[ig3]);
    ind_ran++;
    r = ir1 * NormRANu;
    // printf("r=%f\n",r);
    return r;
}
///Fin números aleatorios (eso si, acordarse de inicializar la semilla)










///     --------MAIN--------       ///

int main()
{
    //Variables que se usan mucho en el main por unas razones o por otras
    int V=L*L*L, i, energia;
    double beta=0.72;



    //Declaramos direccionarios
    int xp[L], yp[L], zp[L], xn[L+1], yn[L], zn[L];

    //Inicializamos los direccionarios
    direccionarios(xp, yp, zp, xn, yn, zn);



    /*    Array de spines, como a cada vértice le corresponden d=3 spines, el primer indice identifica el vértice
    y el segundo indice identifica en que dirección queda el spin correspondiente. Siendo x:=0, y:=1, z:=2    */
    int **spines;
    spines = malloc(V * sizeof(int *));
    for ( i = 0; i < V; i++)
        spines[i] = malloc(3 * sizeof(int));  // 3 direcciones para cada vértice




    // Inicializamos la rueda de números random
    ini_ran(123456789);



    //Le damos una configuración inicial, el flag: 0 (random), 1 o -1 (todos iguales con el respectivo valor 1 o -1)
    //Recomendable usar una de todos iguales para temperaturas bajas (betas altas > 0.7613) y random para altas (betas bajas < 0.7613)
    //Aunque, si se quiere ver la transición de fase mejor la random
    int flag;
    flag = 1;
    //printf("Le damos una configuración inicial, el flag: 0 (random), 1 o -1 (todos iguales con el respectivo valor 1 o -1) \nRecomendable usar una de todos iguales para temperaturas bajas (betas altas) y random para altas (betas bajas) \nAunque, si se quiere ver la transición de fase mejor la random");
    //scanf("flag = %d", &flag);
    configuracion_inicial (flag, spines);



    //Calculamos la energía total de la configuración inicial (sin multiplicar por beta)
    energia = energia_total (xp, yp, zp, spines);
    //energia = beta*energia;
    printf("EnergIa inicial: %f\n", energia*beta);



    //COMPROVACIONES (eliminar o comentar en el futuro)
    //Parece que rula bien :) con 2000 transformaciones de gauge no cambia la energía
    /*
    int indice_random;
    for (i = 0; i < 2000; i++){
        indice_random = (int)(V*Random());
        gauge (xn, yn, zn, spines, indice_random);
        energia = energia_total (xp, yp, zp, spines, beta);
        printf("%f\n", energia);
    }
    */



    //Calculamos las probabilidades del cambio
    double probabilidades[5];
    probabilidades[0] = exp(-8*beta);  //Una vez son mayores que uno se acepta siempre el valor así que pela un poco
    probabilidades[1] = exp(-4*beta);
    probabilidades[2] = 1;
    probabilidades[3] = exp(4*beta);
    probabilidades[4] = exp(8*beta);

    /*  Para el caso lim de beta-->infinito (T=0) basta con poner probabilidades[0] = probabilidades[1] = 0
    probabilidades[3] = probabilidades[4] = 1; probabilidades[2] puede ser cualquier valor, pero 0 hace que se llegue
    más rápido al equilibrio.    */



    ///Metropolis, buscar el equilibrio
    //Ejecutamos Metrópolis, recorriendo M veces toda la malla. Aqui el numero de pasos debe ser suficiente para que el
    // sistema termalice. Una vez alcanzado el equilibrio pasamos al tema de los loops de wilson y tal
    int M = 50; // M = NUMERO DE PASOS METROPOLIS
    for (i = 0; i < M; i++){
        Metropolis (xp, yp, zp, xn, yn, zn, spines, probabilidades, &energia);

        //COMPROBACION, ver que en general la energía del sistema va bajando hasta llegar a un equilibrio
        //printf("%d\n", energia);
    }

    //COMPROBACION de que la energia dada por metropolis y la del sistema sean la misma
    /*
    energia = energia_total (xp, yp, zp, spines);
    printf("Comprobacion: %d\n", energia);
    */

    printf("Energia final: %f\n", energia*beta);




    ///LOOPS DE WILSON
    /*
    Una vez alcanzado el equilibrio, se puede empezar a medir los loops de Wilson. Parece que va a haber que hacer muchas
    repeticiones y este va a ser el cuello de botella del trabajo. Automatizar el proceso me parece muy complicado así que,
    es fácil que lo unico que podamos hacer sea representar la energía con respecto de las iteraciones de montecarlo hasta
    que la energía se mantenga cerca de una zona de equilibrio y por tanto se termalice el sistema. Guardamos la
    configuración que salga y hacemos la estadistica con los loops de wilson.

    Oooooo esa sería la forma rigurosa de hacerlo. Si nuestra ética de trabajo es digamos flexible, podemos hacer una
    chapuza que consiste en: Termalizamos una sola vez, calculamos estadística, hacemos pasar a la configuracion ya
    termalizada por unos cuantos pasos de metrópolis, calculamos estadística, metrópolis, estadística...
    */

    /*El objetivo de esta parte es calcular el promedio de los spines en los loops de wilson en funcion del lado (que luego
    dependen del area pero le sabemos a calcular el area de un cuadrado conocido el lado*/

    //Por no hacer un lio con los índices, el indice vale el lado,  y el vector en promedio_wilson[0] vale siempre 0
    double promedio_wilson[lado_max+1];
    for (i = 0; i < lado_max+1; i++){
        promedio_wilson[i] = 0;
    }

    estadistica_wilson(xp, yp, zp, xn, yn, zn, spines, promedio_wilson);


    ///MÉTODO CHAPUZAS DE HACER ESTADÍSTICA
    int j;
    for (i = 1; i < 10; i++){
        for (j = 1; j < M; j++){
            Metropolis (xp, yp, zp, xn, yn, zn, spines, probabilidades, &energia);
            estadistica_wilson(xp, yp, zp, xn, yn, zn, spines, promedio_wilson);
        }
    }



    for (i = 1; i < lado_max+1; i++){
        printf("Lado = %d: promedio de wilson = %f \n", i, promedio_wilson[i]);
    }

    //Lo pasamos a archivo para ver en gnuplot
    int x;
    double y;
    FILE *Fconfig;
    Fconfig=fopen("promedio_wilson.dat","wt");
    for(i = 1; i < 11; i++){
        //Para linealizar: x = n^2, y = ln(promedio_wilson); --> y = cte - sigma*x
        x = i*i;
        y = log(promedio_wilson[i]);
        fprintf(Fconfig,"%d %f\n", x, y);
    }
    fclose(Fconfig);




    //Liberamos la memoria
    for (i = 0; i < V; i++)
        free(spines[i]);
    free(spines);
}

///     ------FIN DEL MAIN------     ///









//Definimos la forma medio tramposa de movernos por un array 1D como si fuese 3D
void direccionarios(int *xp, int *yp, int *zp, int *xn, int *yn, int *zn){

    int i, L_cuadrado;
    L_cuadrado = L*L;

    for (i = 0; i < (L - 1); i++){
        xp[i] = 1;
        yp[i] = L;
        zp[i]= L_cuadrado;
        xn[i+1] = -1;
        yn[i+1] = -L;
        zn[i+1]= -L_cuadrado;
    }

    xp[L - 1] = -(L - 1);
    yp[L - 1] = -L * (L - 1);
    zp[L - 1] = -L_cuadrado * (L - 1);
    xn[0] = (L - 1);
    yn[0] = L * (L - 1);
    zn[0] = L_cuadrado * (L - 1);
}




/*
A cada vértice le corresponden 3 plaquetas, tantas como planos perpendiculares entre si se puedan definir;
se calcula como combinaciones de d elementos tomadas de 2 en 2.
Por tanto, 3 funciones para las 3 plaquetas correspondientes a cada plano
*/
///OGT: Aunque las llame energías, en realidad no lo son pues NO están multiplicadas por beta
int energia_plaqueta_plano_xy (int *xp, int *yp, int **spines, int n, int i, int j){
    int producto;
    producto = -spines[n][0]*spines[n][1]*spines[n+xp[i]][1]*spines[n+yp[j]][0];
    return producto;
}

int energia_plaqueta_plano_xz (int *xp, int *zp, int **spines, int n, int i, int k){
    int producto;
    producto = -spines[n][0]*spines[n][2]*spines[n+xp[i]][2]*spines[n+zp[k]][0];
    return producto;
}

int energia_plaqueta_plano_yz (int *yp, int *zp, int **spines, int n, int j, int k){
    int producto;
    producto = -spines[n][1]*spines[n][2]*spines[n+yp[j]][2]*spines[n+zp[k]][1];
    return producto;
}


//Recorremos toda la red sumando la energía de cada plaqueta, y ahora si multiplicamos por beta
int energia_total (int *xp, int *yp, int *zp, int **spines){
    unsigned int i, j, k, n;
    int E;

    n=0;
    E=0;
    for (k = 0; k < L; k++){
        for (j = 0; j < L; j++){
            for (i = 0; i < L; i++){
                E += energia_plaqueta_plano_xy (xp, yp, spines, n, i, j);
                E += energia_plaqueta_plano_xz (xp, zp, spines, n, i, k);
                E += energia_plaqueta_plano_yz (yp, zp, spines, n, j, k);
                n++;
            }
        }
    }


    return E;
}



//El flag marca si queremos una configuracion con spines random o uniforme
void configuracion_inicial (int flag, int **spines){
    int V=L*L*L, i, j;

    //Todos 1
    if (flag == 1){
        for (i = 0; i < V; i++){
            for (j = 0; j < 3; j++){
                spines[i][j]=1;
            }
        }
    }

    //Todos -1
    if (flag == -1){
        for (i = 0; i < V; i++){
            for (j = 0; j < 3; j++){
                spines[i][j]=-1;
            }
        }
    }

    //Random
    if (flag == 0){
        for (i = 0; i < V; i++){
            for (j = 0; j < 3; j++){
                if (Random() < 0.5)
                    spines[i][j] = 1;
                else
                    spines[i][j] = -1;
            }
        }
    }

    //Esto no es muy importante, pero genera un fichero de la configuración inicial
    FILE *Fconfig;
    Fconfig=fopen("conf_inicial.dat","wt");
    for(i=0;i<V;i++){
        for (j = 0; j < 3; j++)
            fprintf(Fconfig,"%d ",spines[i][j]);
    }
    fclose(Fconfig);
}



/*
COMPROBACION
Esta función hace una transformada de gauge, su función solo es comprobar que todo rule bien, se puede
eliminar en un futuro
Si rula bien, independientemente de que n elijamos, la energía no debería cambiar
*/
void gauge (int *xn, int *yn, int *zn, int **spines, int n){
    int i, mx, my, mz, aux;

    mx = n%L;
    aux= (n-mx)/L;
    my = aux%L;
    mz = ((aux-my)/L)%L;


    //COMPROBACION de la comprobación, xd esto es muy meta
    /*
    printf("Antes de gauge: ");
    for (i=0;i<3;i++){
        printf("%d ", spines[n][i]);
    }
    printf("%d ", spines[n+xn[mx]][0]);
    printf("%d ", spines[n+yn[my]][1]);
    printf("%d ", spines[n+zn[mz]][2]);
    printf("\n");


    printf("%d ", mx);
    printf("%d ", my);
    printf("%d ", mz);
    printf("\n");
    */


    for (i=0;i<3;i++){
        spines[n][i] = -spines[n][i];
    }
    spines[n+xn[mx]][0] = -spines[n+xn[mx]][0];
    spines[n+yn[my]][1] = -spines[n+yn[my]][1];
    spines[n+zn[mz]][2] = -spines[n+zn[mz]][2];


    //COMPROBACION de la comprobación xd
    /*
    printf("Después de gauge: ");
    for (i=0;i<3;i++){
        printf("%d ", spines[n][i]);
    }
    printf("%d ", spines[n+xn[mx]][0]);
    printf("%d ", spines[n+yn[my]][1]);
    printf("%d ", spines[n+zn[mz]][2]);
    printf("\n");
    */
}



///Final de sacar solo la parte de "Ising normal", un metropolis para llegar al equilibrio térmico
void Metropolis (int *xp, int *yp, int *zp, int *xn, int *yn, int *zn, int **spines, double *probabilidades, int *energia){
    unsigned int i, j, k, n;
    int contribucion_energia, Indice;


    n = 0;
    for (k = 0; k < L; k++){
        for (j = 0; j < L; j++){
            for (i = 0; i < L; i++){

                /*
                Esto es un puto lio, la bromita del Tarancon de querer que todo vaya en un array 1D hace que haya
                que ir cambiando los indices (i,j,k) cuando llamamos a las energia_plaqueta_[]. Como solo cambian
                en una unidad (a menos que estemos en el borde) el vector de cambio es xn (o xp).
                */

                //COMPROBACION: De tanto lio que era tocaba hacer muchas comprobaciones xd
                /*
                printf("Indices: %d %d %d\n", n, n+yn[j], n+zn[k]);
                //int compro=n+zn[k];
                //printf("%d \n", spines[compro+zp[k_primo]][0]);
                printf("%d ", energia_plaqueta_plano_xy(xp, yp, spines, n, i, j));
                printf("%d ", energia_plaqueta_plano_xy(xp, yp, spines, n+yn[j], i, j+xn[j]));
                printf("%d ", energia_plaqueta_plano_xz(xp, zp, spines, n, i, k));
                printf("%d\n", energia_plaqueta_plano_xz(xp, zp, spines, n+zn[k], i, k+xn[k]));
                */


                ///Lo importante de la función
                //Calculamos la energía correspondiente a un spin, es decir, a las 4 plaquetas que toca
                //Un cambio por cada spin, y  a cada vértice le corresponden 3
                contribucion_energia = energia_plaqueta_plano_xy(xp, yp, spines, n, i, j)+energia_plaqueta_plano_xy(xp, yp, spines, n+yn[j], i, j+xn[j])+energia_plaqueta_plano_xz(xp, zp, spines, n, i, k)+energia_plaqueta_plano_xz(xp, zp, spines, n+zn[k], i, k+xn[k]);
                //printf("\n%d %d %d %d", n, n+yn[j], n+zn[k], Indice);
                Indice = (int)(contribucion_energia/2+2);
                if (Random() < probabilidades[Indice]){
                    spines[n][0] = -spines[n][0];
                    *energia += -2*contribucion_energia;
                }
                //printf("E = %d\n", *energia);

                contribucion_energia = energia_plaqueta_plano_xy(xp, yp, spines, n, i, j)+energia_plaqueta_plano_xy(xp, yp, spines, n+xn[i], i+xn[i], j)+energia_plaqueta_plano_yz(yp, zp, spines, n, j, k)+energia_plaqueta_plano_yz(yp, zp, spines, n+zn[k], j, k+xn[k]);
                Indice = (int)(contribucion_energia/2+2);
                if (Random() < probabilidades[Indice]){
                    spines[n][1] = -spines[n][1];
                    *energia += -2*contribucion_energia;
                }
                //printf("E = %d\n", *energia);

                contribucion_energia = energia_plaqueta_plano_xz(xp, zp, spines, n, i, k)+energia_plaqueta_plano_xz(xp, zp, spines, n+xn[i], i+xn[i], k)+energia_plaqueta_plano_yz(yp, zp, spines, n, j, k)+energia_plaqueta_plano_yz(yp, zp, spines, n+yn[j], j+xn[j], k);
                Indice = (int)(contribucion_energia/2+2);
                if (Random() < probabilidades[Indice]){
                    spines[n][2] = -spines[n][2];
                    *energia += -2*contribucion_energia;
                }
                //printf("E = %d\n", *energia);


                n++;
            }
        }
    }
}



///Funciones Loops de Willson
//Recorremos un determinado loop cerrado por orden.
//Tiene pinta de que no hacia falta 3 funciones distintas, pero yo creo que así es menos lio y probablemente más eficiente
int loop_wilson_xy(int *xp, int *yp, int *xn, int *yn, int **spines, int n, int i, int j, int lado){
    int wilson = 1, t;

    for(t = 0; t < lado; t++){
        wilson = wilson*spines[n][1];
        n += yp[j];
        j += xp[j];
    }

    for(t = 0; t < lado; t++){
        wilson = wilson*spines[n][0];
        n += xp[i];
        i += xp[i];
    }

    for(t = 0; t < lado; t++){
        n += yn[j];
        wilson = wilson*spines[n][1];
        j += xn[j];
    }

    for(t = 0; t < lado; t++){
        n += xn[i];
        wilson = wilson*spines[n][0];
        i += xn[i];
    }

    return wilson;
}

int loop_wilson_xz(int *xp, int *zp, int *xn, int *zn, int **spines, int n, int i, int k, int lado){
    int wilson = 1, t;

    for(t = 0; t < lado; t++){
        wilson = wilson*spines[n][2];
        n += zp[k];
        k += xp[k];
    }

    for(t = 0; t < lado; t++){
        wilson = wilson*spines[n][0];
        n += xp[i];
        i += xp[i];
    }

    for(t = 0; t < lado; t++){
        n += zn[k];
        wilson = wilson*spines[n][2];
        k += xn[k];
    }

    for(t = 0; t < lado; t++){
        n += xn[i];
        wilson = wilson*spines[n][0];
        i += xn[i];
    }

    return wilson;
}

//El problema de hacerlo en una sola funcion es esta, que depende también de xp, xn.
//Seria solo poner un nuevo vector pero bueno, ya da palo
int loop_wilson_yz(int *xp, int *yp, int *zp, int *xn, int *yn, int *zn, int **spines, int n, int j, int k, int lado){
    int wilson = 1, t;

    for(t = 0; t < lado; t++){
        wilson = wilson*spines[n][1];
        n += yp[j];
        j += xp[j];
    }

    for(t = 0; t < lado; t++){
        wilson = wilson*spines[n][2];
        n += zp[k];
        k += xp[k];
    }

    for(t = 0; t < lado; t++){
        n += yn[j];
        wilson = wilson*spines[n][1];
        j += xn[j];
    }

    for(t = 0; t < lado; t++){
        n += zn[k];
        wilson = wilson*spines[n][2];
        k += xn[k];
    }

    return wilson;
}



//Calculamos el valor esperado de los loops de
void estadistica_wilson(int *xp, int *yp, int *zp, int *xn, int *yn, int *zn, int **spines, double *promedio_wilson){
    unsigned int lado, i, j, k, n, V = L*L*L;

    /*    No estoy muy seguro de porque, pero dijo en clase que el lado máximo de un loop de wilson era la mitad de L por las
    condiciones toroidales. Según mis "calculos" valdria L-1, pero vamos a hacerle caso a follana que me parece a mi que
    sabe más cosas que nosotros xd     */


    /*    hay que hacer mucha estadística, así que para cada configuracion calculamos el promedio a todos los loops de
    wilson cuadrados posibles con lado =< lado_max.

    He vuelto a leer los apuntes y dice como de que probemos hasta loops de máximo 10*10, bueno se ignoran los de mas o se
    */

    for (lado = 1; lado < lado_max+1; lado++){
        n = 0;
        for (k = 0; k < L; k++){
            for (j = 0; j < L; j++){
                for (i = 0; i < L; i++){
                    promedio_wilson[lado] += loop_wilson_xy(xp, yp, xn, yn, spines, n, i, j, lado);
                    promedio_wilson[lado] += loop_wilson_xz(xp, zp, xn, zn, spines, n, i, k, lado);
                    promedio_wilson[lado] += loop_wilson_yz(xp, yp, zp, xn, yn, zn, spines, n, j, k, lado);

                    n++;
                }
            }
        }

        promedio_wilson[lado] = promedio_wilson[lado]/V;

    }

}
