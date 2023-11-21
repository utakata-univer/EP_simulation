#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ITERATION 1000000
#define FILENAME "Langevin-1D.txt"
#define DATA 10000

double random_Gauss();

int main(){
		
		double t, dt;
			double x, dx, vx;
				int i;
					double gamma, a;
						FILE* fp;
							int s_rate = ITERATION / DATA;
								
								dt = 0.0001;
									x = 0.0;
										gamma = 1.0;
											a = 1.0;
												t = 0.0;
													
													fp = fopen(FILENAME, "wt");
														
														for(i = 0; i < ITERATION; i++){
																	
																	if(i % s_rate == 0){
																					fprintf(fp,"%lf\t%lf\n", t, x);
																							}
																			vx = (sqrt(2.0 * a * gamma / dt) * random_Gauss()) / gamma;
																					dx = dt * vx;	
																							x = x + dx;
																									t = t + dt;
																										}
															fclose(fp);
}

double random_Gauss(){
		static int flag = 0;
			
			double r;
				double ran1;
					double ran2;
						
						if (flag == 0) {
									srand((unsigned int)time(NULL));
											flag = 1;
												}
							
							ran1 = ((double)(rand())) / ((double)RAND_MAX);
								ran2 = ((double)(rand())) / ((double)RAND_MAX);
									
									r = sqrt(-2.0 * log(ran1)) * cos(2.0 * M_PI * ran2);

										return r;
}



