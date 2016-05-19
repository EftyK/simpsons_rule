/* File:    omp_regra_simpson.c
 * Purpose: Aproximar a integral definida pela área sob arcos de parábola que interpolam a função.
 *
 * Input:   x1, x2, n, a, b, (c)
 * Output:  Estimacao do integral de uma funcao
 *
 * Compile: gcc -g -Wall -fopenmp -o omp_regra_simpson omp_regra_simpson.c -lm
 * Usage:   ./omp_regra_simpson <number of threads>
 *
 * Notes:   
 *   O usuario pode escolher o tipo e as constantes da funcao (funcoes semi-personalizadas)
 *  
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void Usage(char* prog_name);
double f(double x, int num);   /* Function we're integrating */
double Simpson(double a, double b, int n, int thread_count, int num);
int escolher_func(void);

int flag = 0; /* Esta 'flag' e usada para pedir os valores dos coeficientes apenas uma vez */

int main(int argc, char* argv[]) {
   double  global_result = 0.0;  /* Store result in global_result */
   double  x1, x2;                 /* Left and right endpoints      */
   int     n;                    /* Total number of trapezoids    */
   int     thread_count;
   int     chosen_function;
  

   if (argc != 2) Usage(argv[0]);
   thread_count = strtol(argv[1], NULL, 10);
   printf("Entre com os intervalos de integracao x1, x2\n");
   scanf("%lf %lf", &x1, &x2);
   
   do{                   /* protecao caso o usuario digite um valor errado */
   printf("Entre com o numero de subintervalos (multiplo de 2)\n");
   scanf("%d", &n);
   }while((n%2) != 0);
   
   chosen_function=escolher_func();
   
   global_result = Simpson(x1, x2, n, thread_count, chosen_function);

   printf("With n = %d intervals, our estimate\n", n);
   printf("of the integral from %f to %f = %.14e\n",
      x1, x2, global_result);
   return 0;
}  /* main */

/*--------------------------------------------------------------------
 * Function:    Usage
 * Purpose:     Print command line for function and terminate
 * In arg:      prog_name
 */
void Usage(char* prog_name) {

   fprintf(stderr, "usage: %s <number of threads>\n", prog_name);
   exit(0);
}  /* Usage */
/*------------------------------------------------------------------
 * Function:    ecolher_func
 * Purpose:     Compute value of function to be integrated
 * Input arg:   x
 * Return val:  f(x)
 */
int escolher_func(void) {

    int num;
    
do{                   /* protecao caso o usuario digite um valor errado */
   printf("Escolha qual tipo de funcao voce deseja integrar\n");
   printf("1 - f(x) = ax^3+bx^2+cx+d\n2 - t(x) = e^(a*x)*cos(b*x)\n3 - v(x) = a/(x+b)\n4 - k(x) = a*seno(b*x)\n");
   printf("5 - r(x) = a*e^(b*x)\n6 - g(x) = a*log(b*x)\n");
   scanf("%d", &num);
}while(num<1 || num>6);
   
   return num;
} /*escolher func*/
/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input arg:   x
 * Return val:  f(x)
 */
double f(double x, int num){
double return_val;
static double a, b, c, d;

   switch(num){
      case 1:
         if(flag == 0){
	 flag = 1;
	 printf("\nVoce escolheu a funcao f(x) = ax�+bx�+cx+d \n");
         printf("Entre com os valores de 'a', 'b', 'c' e 'd'\n");
         scanf("%lf %lf %lf %lf", &a, &b, &c, &d);
         }
	 return_val = a*x*x*x + b*x*x + c*x + d;
         return return_val;
	 break;

      case 2:
         if(flag == 0){
	 flag = 1;
         printf("\nVoce escolheu a funcao t(x) = e^(a*x)*cos(b*x)\n");
         printf("Entre com os valores de 'a' e 'b'\n");
         scanf("%lf %lf", &a, &b);
	 }
         return_val = (exp(a*x))*cos(b*x);
         return return_val;
      break;

      case 3:
         if(flag == 0){
	 flag = 1;
         printf("\nVoce escolheu a funcao v(x) = a/(x+b)\n");
         printf("Entre com os valores de 'a' e 'b'\n");
         scanf("%lf %lf", &a, &b);
	 }
         return_val = a/(x+b);
         return return_val;
      break;

      case 4:
         if(flag == 0){
	 flag = 1;
         printf("\nVoce escolheu a funcao k(x) = a*seno(b*x)\n");
         printf("Entre com os valores de 'a' e 'b'\n");
         scanf("%lf %lf", &a, &b);
	 }
         return_val = a*sin(b*x);
         return return_val;
      break;
      case 5:
         if(flag == 0){
	 flag = 1;
         printf("\nVoce escolheu a funcao r(x) = a*e^(b*x)\n");
         printf("Entre com os valores de 'a' e 'b'\n");
         scanf("%lf %lf", &a, &b);
	 }
         return_val = a*exp(b*x);
         return return_val;
      break;
      case 6:
         if(flag == 0){
	 flag = 1;
	 printf("\nVoce escolheu a funcao g(x) = a*log(b*x) \n");
         printf("Entre com os valores de 'a' e 'b'\n");
         scanf("%lf %lf", &a, &b);
         }
	 return_val =a*log(b*x);
         return return_val;
      break;
   }
return 0;
} /*f*/

/*------------------------------------------------------------------
 * Function:    Simpson
 * Purpose:     Use Simpson rule to estimate definite integral
 * Input args:  
 *    x1: left endpoint
 *    x2: right endpoint
 *    n: number of trapezoids
 * Return val:
 *    integral:  estimate of integral from a to b of f(x)
 */

double Simpson(double x1, double x2, int n,int thread_count,int num){
 
double h;
double x,y;
int i;
double integral;
double start, elapsed;

h = (x2-x1)/n; 

integral = f(x1,1);

start = omp_get_wtime();
#  pragma omp parallel for num_threads(thread_count) \
      reduction(+: integral) 
for(i=1; i<n; i++){
   x=x1+h*i;
   y=f(x,num);
   if(i%2 == 0)                         /*Resto da divisao por 2 nulo*/
      integral = integral + 2*y;
   else
      integral = integral + 4*y;
   }
elapsed = omp_get_wtime()- start;  
  
integral = (h/3)*(integral + f(x2,num));

printf("Time elapsed = %f\n", elapsed);

return integral;
}/*Simpson*/