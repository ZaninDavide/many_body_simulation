#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#define CROW 4 // Number of cells in a row
#define M 4 // Number of particles in a cell
#define m 1.0
#define Kb 1.0
#define epsilon 1.0
#define PI 3.1415926535897932384626433
#define BINS 150 // number of bins for the histogram of g(r)

#define N (CROW*CROW*CROW*M) // Total number of particles
double T = 0.0;
double rho = 0.0; // fluid density
double L = 0.0; // periodicity, will be overwritten

double lennard_jones_L_halves = 0.0;
double* histogram;


// Single particle
double sp_delta1 = 0.0; // Length of the deterministic jump (proportional to the forces)
double sp_delta2 = 0.1; // Length of the random jump (proportional to a random direction)
// Many particles
double mp_delta1 = 0.0; // Length of the deterministic jump (proportional to the forces)
double mp_delta2 = 0.005; // Length of the random jump (proportional to a random direction)

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
    double ener; // energia
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

void add_to_distribution_histogram(Particle S[]) {
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
            } else {
                histogram[bin_index] += 1; 
            }
        }
    }
}

void save_distribution_histogram(FILE* file, unsigned int steps) {
    double bin_width = L / BINS;
    // Normalize and print to file
    for (int j = 1; j < BINS; j++) {
        double dV = 4*PI/3.0*((j+1)*(j+1)*(j+1) - j*j*j)*bin_width*bin_width*bin_width;
        histogram[j] /= dV * rho * N * (steps/2);
        fprintf(file, "%10.5e %10.5e\n", (j + 0.0) / BINS, histogram[j]);
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

    printf("Step: 0/%d", steps);
    for(int i = 1; i <= steps; i++) {
        if(i % 100 == 0) { printf("\rStep: %d/%d", i, steps); fflush(stdout); }
        // Recenter the positions
        // recenter_particles(S0);
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
        if(histogram && i >= steps/2) {
            add_to_distribution_histogram(S1);
        }
        // Preparo il prossimo ciclo
        swap_pointer((void**)&F0, (void**)&F1);
        swap_pointer((void**)&S0, (void**)&S1);
    }
    printf("\rDone!\n");

    if(steps % 2 == 1) swap_pointer((void**)&S0, (void**)&S1);

    free(F0);
    free(F1);
    free(S1);

    return measures;
}

double lennard_jones_interaction(Particle S[], int i, int j) {
    // Calculate V(L/2)
    vec3 Rij = (vec3) {
        (S[i].x - S[j].x) - L * rint((S[i].x - S[j].x)/L),
        (S[i].y - S[j].y) - L * rint((S[i].y - S[j].y)/L),
        (S[i].z - S[j].z) - L * rint((S[i].z - S[j].z)/L),
    };
    double r2 = Rij.x*Rij.x + Rij.y*Rij.y + Rij.z*Rij.z;
    if(r2 < (L/2)*(L/2)) { // truncate forces after r = L/2
        double sr6 = pow(r2, -3.0);
        return 4*epsilon*sr6*(sr6 - 1) - lennard_jones_L_halves;
    }
    return 0;
}

double get_potential_energy(Particle S[], double (*get_interaction)(Particle[], int, int)) {
    // Calculate total potental energy
    double Epot = 0.0;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < i; j++) {
            Epot += get_interaction(S, i, j);
        }
    }
    return Epot;
}

double move_one_particle(Particle S0[], double (*get_interaction)(Particle[], int, int), int k, vec3 delta) {
    double interaction_before = 0.0;
    double interaction_after = 0.0;
    for(int i = 0; i < N; i++) {
        if(i == k) continue;
        interaction_before += get_interaction(S0, k, i);
    }
    S0[k].x += delta.x;
    S0[k].y += delta.y;
    S0[k].z += delta.z;
    for(int i = 0; i < N; i++) {
        if(i == k) continue;
        interaction_after += get_interaction(S0, k, i);
    }
    return interaction_after - interaction_before;
}

// Get random number uniformly distributed between min and max
double sample_uniform(double min, double max) {
    double x = rand()/(RAND_MAX + 1.0);
    return min + x*(max - min);
}

double V0 = 0.0; // V0 is global to be used by measure()
// Simulate the system with M(RT)^2
Measure* metropolis(Particle S0[], double (*get_interaction)(Particle [], int, int), void (*get_forces)(Particle S[], vec3 F[]), Measure (*measure)(Particle S[], vec3 F[], double t), unsigned int steps, bool single_particle_moves, FILE* file) {
    vec3* F = malloc(N * sizeof(vec3));
    Particle* S1 = NULL;
    if(single_particle_moves == false) { S1 = malloc(N * sizeof(Particle)); }

    Particle* toFree = S1;

    double V1 = 0.0;
    V0 = 0.0;

    // Scrivo in output lo stato iniziale
    write_particles_to_file(S0, 0, file);

    // Inizializzo le forze
    get_forces(S0, F);

    // Inizializzo l'energia potenziale
    V0 = get_potential_energy(S0, get_interaction);

    // Inizializzo le misure e misuro
    Measure* measures = NULL;
    if(measure) {
        measures = calloc(steps + 1, sizeof(Measure));
        measures[0] = measure(S0, F, 0.0);
    }

    int jumps = 0;
    printf("Step: 0/%d", steps);
    for(int i = 1; i <= steps; i++) {
        if(single_particle_moves == true && i % 100 == 0) { 
            printf("\rStep: %d/%d, Acceptance: %.2f", i, steps, jumps / (double)i / N); fflush(stdout); 
        }else if(single_particle_moves == false && i % 100 == 0){ 
            printf("\rStep: %d/%d, Acceptance: %.2f", i, steps, jumps / (double)i); fflush(stdout);
        }
        // Genero una nuova configurazione S1
        if(single_particle_moves == true) {
            // Muovo una particella alla volta, e ad ogni spostamento provo un salto
            // Qui uso solo S0 dato che coinciderebbe sempre con S1
            for(int k = 0; k < N; k++){
                vec3 delta = (vec3){
                    sp_delta1*F[k].x + sp_delta2*sample_uniform(-1,1),
                    sp_delta1*F[k].y + sp_delta2*sample_uniform(-1,1),
                    sp_delta1*F[k].z + sp_delta2*sample_uniform(-1,1),
                };
                double deltaV = move_one_particle(S0, get_interaction, k, delta);
                // Calcolo la probabilità di salto alla nuova configurazione
                V1 = V0 + deltaV;
                double p_ratio = exp(-(V1 - V0)/Kb/T);
                double a = (p_ratio > 1) ? 1 : p_ratio;
                if(sample_uniform(0, 1) < a) {
                    // Accetto il salto quindi preparo il potenziale per il prossimo ciclo
                    V0 = V1;
                    jumps += 1;
                } else {
                    // Inverto il movimento dato che rifiuto il salto
                    S0[k].x -= delta.x;
                    S0[k].y -= delta.y;
                    S0[k].z -= delta.z;
                }
            }
            get_forces(S0, F);
        }else{
            // Muovo tutte le particelle in un colpo solo
            for(int k = 0; k < N; k++){
                S1[k].x = S0[k].x + mp_delta1*F[k].x + mp_delta2*sample_uniform(-1,1);
                S1[k].y = S0[k].y + mp_delta1*F[k].y + mp_delta2*sample_uniform(-1,1);
                S1[k].z = S0[k].z + mp_delta1*F[k].z + mp_delta2*sample_uniform(-1,1);
            }
            // Calcolo la probabilità di salto alla configurazione S1
            V1 = get_potential_energy(S1, get_interaction);
            double p_ratio = exp(-(V1 - V0)/Kb/T);
            double a = (p_ratio > 1) ? 1 : p_ratio;
            // Se devo saltare salto
            if(sample_uniform(0, 1) < a) {
                // Preparo le forze e il potenziale per il prossimo ciclo
                get_forces(S1, F);
                V0 = V1;
                // La posizione per il prossimo ciclo diventa S1 
                swap_pointer((void**)&S0, (void**)&S1);
                jumps += 1;
            } 
        }
        // Ora, indipendentemente da se ho saltato o meno, S0 contiene le prossime posizioni (e velocità anche se in metropolis non hanno senso le v)
        // Scrivo in output le nuove posizioni
        write_particles_to_file(S0, i, file);
        // Misura delle osservabili
        if(measure && measures) {
            measures[i] = measure(S0, F, i);
        }
        if(histogram && i >= steps/2) {
            add_to_distribution_histogram(S0);
        }
    }
    if(single_particle_moves == true){ 
        printf("\rAcceptance: %f\n", jumps / ((double) steps) / N); 
    } else { 
        printf("\rAcceptance: %f\n", jumps / ((double) steps)); 
    }

    free(F);
    free(toFree);

    return measures;
}

void initialize_positions(Particle S0[]) {
    double a = L / (double)CROW;
    
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

        double y1 = rand()/(RAND_MAX + 1.0);
        double y2 = rand()/(RAND_MAX + 1.0);
        double factor_y = sigma*sqrt(-2*log(1-y2));

        double z1 = rand()/(RAND_MAX + 1.0);
        double z2 = rand()/(RAND_MAX + 1.0);
        double factor_z = sigma*sqrt(-2*log(1-z2));
        
        S0[i].vx = factor_x*cos(2*PI*x1);
        S0[i].vy = factor_y*cos(2*PI*y1);
        S0[i].vz = factor_z*cos(2*PI*z1);

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

double W_forces = 0.0;
// Calculates the forces produced by the Lennard-Jones potential
void lennard_jones(Particle S[], vec3 F[]) {
    double sigma2 = 1.0;

    W_forces = 0.0;
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
                double factor = 24.0*epsilon*sr6*(2.0*sr6 - 1.0) / r2; 
                vec3 Fij = (vec3) { Rij.x * factor, Rij.y * factor, Rij.z * factor };
                F[i].x += Fij.x;
                F[i].y += Fij.y;
                F[i].z += Fij.z;
                F[j].x -= Fij.x;
                F[j].y -= Fij.y;
                F[j].z -= Fij.z;
                W_forces += Fij.x*Rij.x + Fij.y*Rij.y + Fij.z*Rij.z;
            }
        }
    }
    W_forces /= N;
}

FILE* file_measures = NULL;
Measure measure(Particle S[], vec3 F[], double time) {
    // Temperatura <Ek> = 3/2.Kb.T
    double Ek = 0.0;
    for(int k = 0; k < N; k++) {
        Ek += 0.5 * m * (S[k].vx*S[k].vx + S[k].vy*S[k].vy + S[k].vz*S[k].vz);
    }
    double temperature = 2.0 / 3.0 / Kb * (Ek / N);

    // Potential energy
    double Epot = get_potential_energy(S, lennard_jones_interaction);
    
    // Compressibilità
    double compressibilità = 1 + W_forces / 3.0 / rho / Kb / temperature;

    if(file_measures) fprintf(file_measures, "%10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n", time, temperature, Ek / N, Epot / N, (Ek + Epot) / N, compressibilità);

    return (Measure) { temperature, compressibilità, (Ek + Epot) / N };
}

Measure measure_metropolis(Particle S[], vec3 F[], double time) {
    // Potential energy
    double Epot = V0; // get_potential_energy(S, lennard_jones_interaction);
    
    // Compressibilità
    double compressibilità = 1 + W_forces / 3.0 / rho / Kb / T;

    if(file_measures) fprintf(file_measures, "%10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n", time, T, 0.0, Epot / N, (0.0 + Epot) / N, compressibilità);

    return (Measure) { T, compressibilità, (0.0 + Epot) / N };
}

Measure averages(Measure* measures, int start, int end) {
    Measure avg = (Measure) {0, 0, 0};
    for(int i = start; i < end; i++) { 
        avg.temp += measures[i].temp; 
        avg.comp += measures[i].comp; 
        avg.ener += measures[i].ener; 
    }
    avg.temp /= end - start;
    avg.comp /= end - start;
    avg.ener /= end - start;
    return avg;
}

Measure variances(Measure* measures, Measure avg, int start, int end) {
    Measure var = (Measure) {0, 0, 0};
    for(int i = start; i < end; i++) { 
        var.temp += (measures[i].temp - avg.temp)*(measures[i].temp - avg.temp); 
        var.comp += (measures[i].comp - avg.comp)*(measures[i].comp - avg.comp); 
        var.ener += (measures[i].ener - avg.ener)*(measures[i].ener - avg.ener); 
    }
    var.temp /= end - start;
    var.comp /= end - start;
    var.ener /= end - start;
    return var;
}

void averages_and_variances(Measure *measures, int steps, FILE* varAvgBFile) {
    // Global Averages (on the second half of the data)
    Measure avg = averages(measures, steps/2, steps);
    Measure var = variances(measures, avg, steps/2, steps); // Meaningful only with vverlet

    printf("Temperatura = %f ± %f\n", avg.temp, sqrt(var.temp));
    printf("Compressibilità = %f ± %f\n", avg.comp, sqrt(var.comp));
    printf("Energia = %f ± %2.5e\n", avg.ener, sqrt(var.ener));

    // Variances for metropolis
    int start = steps / 2;
    int NN = steps - start; // Size of the data
    for(int B = 1; B < NN / 2; B++) {  // Loop over the number of points per block
        int NB = NN / B; // The number of blocks
        Measure* avgB = malloc(sizeof(Measure) * NB); // Averages in the blocks 
        for(int b = 0; b < NB; b++){
            avgB[b] = averages(measures, start + B*b, start + B*(b + 1));
        }
        Measure varAvgB = variances(avgB, averages(avgB, 0, NB), 0, NB); // Variances in the blocks
        if(varAvgBFile) fprintf(varAvgBFile, "%d %f %f\n", B, sqrt(varAvgB.ener / NB), sqrt(varAvgB.comp / NB));
        free(avgB);
    }

}

void initialize_constants() {
    L = CROW * pow(M / rho, 1.0/3.0); // rho = N/L3 = P3*M/L3

    // Calculate V(L/2)
    double sigma2 = 1.0;
    double sigma_su_L_mezzi_2 = sigma2 / (L/2) / (L/2);
    double sigma_su_L_mezzi_6 = sigma_su_L_mezzi_2 * sigma_su_L_mezzi_2 * sigma_su_L_mezzi_2;
    lennard_jones_L_halves = 4*epsilon*sigma_su_L_mezzi_6*(sigma_su_L_mezzi_6 - 1);
    
    V0 = 0;
    W_forces = 0;
    
    printf("CROW = %d, N = %d, L = %f, T0 = %f, rho = %f\n", CROW, N, L, T, rho);
}

int main() {
    // To use these you need to #define N 1
    // FILE* file_osc = fopen("file_osc.dat", "w+");
    // Particle S0osc[] = { (Particle) { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };
    // vverlet(S0osc, oscillatore_armonico, NULL, 0.01, 5000, file_osc);
    // return 0;

    double simulations[14][3] = {
        // rho, T, sp_delta2
        {1.200, 1.1, 0.10}, 
        {0.010, 1.1, 0.50},
        {0.800, 1.1, 0.13},
        {0.700, 1.1, 0.13},
        {0.600, 1.1, 0.13},
        {0.100, 1.1, 0.50},
        {1.000, 1.1, 0.10},
        {0.750, 1.1, 0.13},
        {0.200, 1.1, 0.50},
        {0.400, 1.1, 0.13},
        {0.001, 1.1, 0.13},
        {0.500, 1.1, 0.13},
        {0.300, 1.1, 0.13},
        {0.250, 1.1, 0.13},
    };

    Particle* S0 = malloc(N * sizeof(Particle));
    histogram = malloc(BINS * sizeof(double));
    for(int k = 0; k < 14; k++){
        char file_measures_name[20]; sprintf(file_measures_name, "measures_%02d.dat", k);
        file_measures = fopen(file_measures_name, "w+");
        memset(histogram, 0, BINS * sizeof(double));

        rho = simulations[k][0]; // Density
        T = simulations[k][1]; // Random Jump Size
        sp_delta2 = simulations[k][2]; // Random Jump Size
        initialize_constants();
        initialize_positions(S0);
        initialize_velocities(S0);

        unsigned int steps = 10000;
        // Measure* measures = vverlet(S0, lennard_jones, measure, 0.001, steps, NULL);
        Measure* measures = metropolis(S0, lennard_jones_interaction, lennard_jones, measure_metropolis, steps, true, NULL);

        char var_avg_B_file_name[20]; sprintf(var_avg_B_file_name, "varAvgB_%02d.dat", k);
        FILE* var_avg_B_file = fopen(var_avg_B_file_name, "w+");
        averages_and_variances(measures, steps, var_avg_B_file);

        // Finalize and export histogram
        // FILE* histogram_file = fopen("histogram.dat", "w+");
        // save_distribution_histogram(histogram_file, steps);

        free(measures);
        fclose(file_measures);
        fclose(var_avg_B_file);
        // fclose(histogram_file);
    }

    // Scatter plot of particles' positions
    /*
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
    */

    return 0;
}