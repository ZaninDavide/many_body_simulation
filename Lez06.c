#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define CROW 2 // Number of cells in a row
#define M 4 // Number of particles in a cell
#define rho 0.8
#define m 1.0
#define Kb 1.0
#define T 1.8
#define epsilon 1.0
#define PI 3.1415926535897932384626433
#define BINS 150 // number of bins for the histogram of g(r)

#define N (CROW*CROW*CROW*M) // Total number of particles
double L = 0.0; // periodicity, will be overwritten

typedef struct Particle {
    double x; double y; double z;
    double vx; double vy; double vz;
} Particle;

typedef struct vec3 {
    double x; double y; double z;
} vec3;

typedef struct Measure {
    double temp; // temperature
    double comp; // compressibilità
} Measure;

void swap_pointer(void** p1, void** p2) {
    void* p3 = *p1;
    *p1 = *p2;
    *p2 = p3;
}

void write_particles_to_file(Particle S[], double t, FILE* file) {
    if(!file) return;
    fprintf(file, "%10.5e ", t);
    for(int i = 0; i < N; i++) {
        fprintf(file, "%10.5e %10.5e %10.5e ", S[i].x, S[i].y, S[i].z);
    }
    fprintf(file, "\n");
}

void recenter_particles(Particle S[]) {
    for (int i = 0; i < N; i++) {
        vec3 pos = (vec3) {
            S[i].x - L * rint(S[i].x/L),
            S[i].y - L * rint(S[i].y/L),
            S[i].z - L * rint(S[i].z/L),
        };
        S[i].x = pos.x;
        S[i].y = pos.y;
        S[i].z = pos.z;
    }
}

// Performs the velocity-Verlet algorithm.
// Returns the list of measurements made at every step of the simulation. 
// Overrides the content of S0.
// Every step of the simulation, the entire state of the system is written to the specified file according to write_particles_to_file
Measure* vverlet(Particle S0[], void (*get_forces)(Particle S[], vec3 F[]), Measure (*measure)(Particle S[], vec3 F[], double t), double dt, unsigned int steps, FILE* file) {
    vec3* F0 = malloc(N * sizeof(vec3));
    vec3* F1 = malloc(N * sizeof(vec3));
    Particle* S1 = malloc(N * sizeof(Particle));

    // Scrivo in output lo stato iniziale
    write_particles_to_file(S0, 0, file);

    // Inizializzo le forze
    get_forces(S0, F0);

    // Inizializzo le misure e misuro
    Measure* measures = NULL;
    if(measure) {
        measures = calloc(steps + 1, sizeof(Measure));
        measures[0] = measure(S0, F0, 0.0);
    }

    for(int i = 1; i <= steps; i++) {
        // Every once in a while recenter the positions
        // if(i % 100 == 100 - 1) { recenter_particles(S0); }
        // Avendo le posizioni, le velocità e le forze al tempo t0 ottengo le nuove posizioni
        for(int k = 0; k < N; k++) {
            S1[k].x = S0[k].x + dt*S0[k].vx + 0.5*F0[k].x/m*dt*dt;
            S1[k].y = S0[k].y + dt*S0[k].vy + 0.5*F0[k].y/m*dt*dt;
            S1[k].z = S0[k].z + dt*S0[k].vz + 0.5*F0[k].z/m*dt*dt;
        }
        // Avendo le nuove posizioni calcolo le nuove forze (qui le nuove velocità non servono)
        get_forces(S1, F1);
        // Avendo le velocità, le forze e le nuove forze posso calcolare le nuove velocità
        for(int k = 0; k < N; k++) {
            S1[k].vx = S0[k].vx + (F0[k].x + F1[k].x)/m/2.0*dt;
            S1[k].vy = S0[k].vy + (F0[k].y + F1[k].y)/m/2.0*dt;
            S1[k].vz = S0[k].vz + (F0[k].z + F1[k].z)/m/2.0*dt;
        }
        // Scrivo in output le nuove posizioni
        write_particles_to_file(S1, i*dt, file);
        // Misura delle osservabili
        if(measure && measures) {
            measures[i] = measure(S1, F1, i*dt);
        }
        // Preparo il prossimo ciclo
        swap_pointer((void**)&F0, (void**)&F1);
        swap_pointer((void**)&S0, (void**)&S1);
    }

    if(steps % 2 == 1) swap_pointer((void**)&S0, (void**)&S1);

    free(F0);
    free(F1);
    free(S1);

    return measures;
}

void initialize_positions(Particle S0[]) {
    double a = L / (double)CROW;

    printf("CROW = %d, N = %d, L = %f, T0 = %f, rho = %f\n", CROW, N, L, T, rho);
    
    vec3 CC[1]  = { (vec3){0, 0, 0} };
    vec3 BCC[2] = { (vec3){0, 0, 0}, (vec3){0.5, 0.5, 0.5} };
    vec3 FCC[4] = { (vec3){0, 0, 0}, (vec3){0.5, 0.5, 0}, (vec3){0.5, 0, 0.5}, (vec3){0, 0.5, 0.5} };
    vec3* CELL_TYPES[] = {NULL, CC, BCC, NULL, FCC};

    unsigned int index = 0;
    for(int nx = 0; nx <= CROW - 1; nx++) {
    for(int ny = 0; ny <= CROW - 1; ny++) {
    for(int nz = 0; nz <= CROW - 1; nz++) {
        for(int p = 0; p < M; p++) {
            S0[index].x = a * nx + a*CELL_TYPES[M][p].x; 
            S0[index].y = a * ny + a*CELL_TYPES[M][p].y; 
            S0[index].z = a * nz + a*CELL_TYPES[M][p].z;
            index += 1;            
        }
    }}}
}

void initialize_velocities(Particle S0[]) {
    srand(9);
    double sigma = sqrt(Kb*T/m);
    for(int i = 0; i < N; i+=2) {
        double x1 = rand()/(RAND_MAX + 1.0);
        double x2 = rand()/(RAND_MAX + 1.0);
        double factor_x = sigma*sqrt(-2*log(1-x2));
        S0[i  ].vx = factor_x*cos(2*PI*x1);

        double y1 = rand()/(RAND_MAX + 1.0);
        double y2 = rand()/(RAND_MAX + 1.0);
        double factor_y = sigma*sqrt(-2*log(1-y2));
        S0[i  ].vy = factor_y*cos(2*PI*y1);

        double z1 = rand()/(RAND_MAX + 1.0);
        double z2 = rand()/(RAND_MAX + 1.0);
        double factor_z = sigma*sqrt(-2*log(1-z2));
        S0[i  ].vz = factor_z*cos(2*PI*z1);

        if(i + 1 < N) {
            S0[i+1].vx = factor_x*sin(2*PI*x1);
            S0[i+1].vy = factor_y*sin(2*PI*y1);
            S0[i+1].vz = factor_z*sin(2*PI*z1);
        }
    }
}

void oscillatore_armonico(Particle S[], vec3 F[]) {
    F[0].x = -S[0].x;
}

// Calculates the forces produced by the Lennard-Jones potential
void lennard_jones(Particle S[], vec3 F[]) {
    double sigma2 = Kb*T/m;
    
    memset(F, 0, N*sizeof(vec3)); // zero out forces
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < i; j++) {
            vec3 Rij = (vec3) {
                (S[i].x - S[j].x) - L * rint((S[i].x - S[j].x)/L),
                (S[i].y - S[j].y) - L * rint((S[i].y - S[j].y)/L),
                (S[i].z - S[j].z) - L * rint((S[i].z - S[j].z)/L),
            };
            double r2 = Rij.x*Rij.x + Rij.y*Rij.y + Rij.z*Rij.z;
            if(r2 < (L/2.0)*(L/2.0)) { // truncate forces after r = L/2
                double sr2 = sigma2 / r2;
                double sr6 = sr2 * sr2 * sr2;
                double factor = 24.0*epsilon*sr6*(2.0*sr6 - 1.0);
                vec3 Fij = (vec3) { Rij.x * factor, Rij.y * factor, Rij.z * factor };
                F[i].x += Fij.x;
                F[i].y += Fij.y;
                F[i].z += Fij.z;
                F[j].x -= Fij.x;
                F[j].y -= Fij.y;
                F[j].z -= Fij.z;
            }
        }
    }
}

FILE* file_measures = NULL;
Measure measure(Particle S[], vec3 F[], double time) {
    double Ek = 0.0;
    for(int k = 0; k < N; k++) {
        Ek += 0.5 * m * (S[k].vx*S[k].vx + S[k].vy*S[k].vy + S[k].vz*S[k].vz);
    }
    // <Ek> = 3/2.Kb.T
    double temperature = 2.0 / 3.0 / Kb * (Ek / N);

    // V(L/2)
    double sigma2 = Kb*T/m;
    double sigma_su_L_mezzi_2 = sigma2 / (L/2) / (L/2);
    double sigma_su_L_mezzi_6 = sigma_su_L_mezzi_2 * sigma_su_L_mezzi_2 * sigma_su_L_mezzi_2;
    double V_L_mezzi = 4*epsilon*sigma_su_L_mezzi_6*(sigma_su_L_mezzi_6 - 1);
    // Potential energy
    double Epot = 0.0;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < i; j++) {
            vec3 Rij = (vec3) {
                (S[i].x - S[j].x) - L * rint((S[i].x - S[j].x)/L),
                (S[i].y - S[j].y) - L * rint((S[i].y - S[j].y)/L),
                (S[i].z - S[j].z) - L * rint((S[i].z - S[j].z)/L),
            };
            double r2 = Rij.x*Rij.x + Rij.y*Rij.y + Rij.z*Rij.z;
            if(r2 < (L/2)*(L/2)) { // truncate forces after r = L/2
                double sr6 = powf(sigma2 / r2, 3.0);
                Epot += 4*epsilon*sr6*(sr6 - 1) - V_L_mezzi;
            }
        }
    }

    // Pressure
    double W = 0.0;
    for(int i = 0; i < N; i++) {
        vec3 Ri = (vec3) {
            S[i].x - L * rint(S[i].x / L),
            S[i].y - L * rint(S[i].y / L),
            S[i].z - L * rint(S[i].z / L),
        };
        W += F[i].x*Ri.x + F[i].y*Ri.y + F[i].z*Ri.z;
    }
    W = W / N;
    double compressibilità = 1 + W / 3.0 / rho / Kb / temperature;

    if(file_measures) fprintf(file_measures, "%10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n", time, temperature, Ek / N, Epot / N, (Ek + Epot) / N, compressibilità);

    return (Measure) { temperature, compressibilità };
}

void distribution_histogram(Particle S[], FILE* file) {
    double histogram[BINS];
    memset(histogram, 0, BINS*sizeof(double)); // zero out histogram
    double bin_width = L / BINS;
    for (int i = 1; i < N; i++) {
        for (int j = 1; j < N; j++) {
            if(i == j) continue;
            vec3 Rij = (vec3) {
                (S[i].x - S[j].x) - L * rint((S[i].x - S[j].x)/L),
                (S[i].y - S[j].y) - L * rint((S[i].y - S[j].y)/L),
                (S[i].z - S[j].z) - L * rint((S[i].z - S[j].z)/L),
            };
            double r = sqrt(Rij.x*Rij.x + Rij.y*Rij.y + Rij.z*Rij.z);
            int bin_index = r / bin_width;
            if(r > sqrt(3)*L) { printf("Distanza fuori scala %f > %f\n", r, L);  }
            if(bin_index > BINS - 1) { 
                printf("Indice fuori scala %d > %d\n", bin_index, BINS); 
                // printf("Si  = [%f, %f, %f]\n", S[i].x, S[i].y, S[i].z); 
                // printf("Sj  = [%f, %f, %f]\n", S[j].x, S[j].y, S[j].z); 
                // printf("Rij = [%f, %f, %f]\n", Rij.x, Rij.y, Rij.z); 
            } else {
                histogram[bin_index] += 1; 
            }
        }
    }

    // normalize and print to file
    for (int j = 1; j < BINS; j++) {
        double dV = 4*PI/3.0*((j+1)*(j+1)*(j+1) - j*j*j)*bin_width*bin_width*bin_width;
        histogram[j] /= dV * rho * N;
        fprintf(file, "%10.5e %10.5e\n", (j + 0.0) / BINS, histogram[j]);
    }
}

void tests() {
    double a_value = 5;
    double b_value = 10;
    double* a = &a_value;
    double* b = &b_value;
    swap_pointer((void**)&a, (void**)&b);
    printf("*a = %f\n", *a);
    printf("*b = %f\n", *b);

    double arr[] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    memset(arr, 0, 5*sizeof(double));
    printf("arr = {%f, %f, %f, %f, %f}\n", arr[0], arr[1], arr[2], arr[3], arr[4]);

    vec3 CC[1]  = { (vec3){0, 0, 0} };
    vec3 BCC[2] = { (vec3){0, 0, 0}, (vec3){0.5, 0.5, 0.5} };
    vec3 FCC[4] = { (vec3){0, 0, 0}, (vec3){0.5, 0.5, 0}, (vec3){0.5, 0, 0.5}, (vec3){0, 0.5, 0.5} };
    vec3* CELL_TYPES[] = {NULL, CC, BCC, NULL, FCC};
    printf("FCC = {  (vec3){%f, %f, %f}, (vec3){%f, %f, %f}, (vec3){%f, %f, %f}, (vec3){%f, %f, %f}  }",
        CELL_TYPES[4][0].x, CELL_TYPES[4][0].y, CELL_TYPES[4][0].z,
        CELL_TYPES[4][1].x, CELL_TYPES[4][1].y, CELL_TYPES[4][1].z,
        CELL_TYPES[4][2].x, CELL_TYPES[4][2].y, CELL_TYPES[4][2].z,
        CELL_TYPES[4][3].x, CELL_TYPES[4][3].y, CELL_TYPES[4][3].z
    );
}

int main() {
    // To use these you need to #define N 1
    // FILE* file_osc = fopen("file_osc.dat", "w+");
    // Particle S0osc[] = { (Particle) { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };
    // vverlet(S0osc, oscillatore_armonico, NULL, 0.01, 5000, file_osc);
    // return 0;

    // tests();
    // return 0;

    file_measures = fopen("measures.dat", "w+");

    // rho = N/L3 = P3*M/L3
    L = CROW * powf(M / rho, 1.0/3.0);

    Particle* S0 = malloc(sizeof(Particle) * N);
    initialize_positions(S0);
    initialize_velocities(S0);

    FILE* file = NULL; // fopen("gas.dat", "w+");
    unsigned int steps = 3000;
    Measure* measures = vverlet(S0, lennard_jones, measure, 0.001, steps, file);

    // Average temperature
    double avg_temp = 0.0;
    double avg_comp = 0.0;
    for(int i = steps / 2; i <= steps; i++) { 
        avg_temp += measures[i].temp; 
        avg_comp += measures[i].comp; 
    }
    avg_temp /= (steps / 2);
    avg_comp /= (steps / 2);
    double sigma_temp = 0.0;
    double sigma_comp = 0.0;
    for(int i = steps / 2; i <= steps; i++) { 
        sigma_temp += (measures[i].temp - avg_temp)*(measures[i].temp - avg_temp); 
        sigma_comp += (measures[i].comp - avg_comp)*(measures[i].comp - avg_comp); 
    }
    sigma_temp /= (steps / 2);
    sigma_comp /= (steps / 2);

    // NON E' NUMERICAMENTE STABILE
    printf("Temperatura = %f ± %f\nCompressibilità = %f ± %f\n", avg_temp, sigma_temp, avg_comp, sigma_comp);

    FILE* histogram = fopen("histogram.dat", "w+");
    distribution_histogram(S0, histogram); // g(r)


    FILE* scatter = fopen("scatter.dat", "w+");
    vec3 center = (vec3) {L/2, L/2, L/2};
    for(int i = 0; i < N; i++) {
        vec3 Ri = (vec3) {
            (S0[i].x - center.x) - L * rint((S0[i].x - L/2)/L),
            (S0[i].y - center.y) - L * rint((S0[i].y - L/2)/L),
            (S0[i].z - center.z) - L * rint((S0[i].z - L/2)/L),
        };
        fprintf(scatter, "%f %f %f\n", Ri.x, Ri.y, Ri.z);
    }

    free(measures);

    return 0;
}