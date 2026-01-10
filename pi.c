/* * Programme de calcul des décimales de PI via l'algorithme de Chudnovsky.
 * Utilise la bibliothèque GMP (assez précis) et Pthreads pour le multithread ainsi que Binary Splitting.
 * Options : 
 * -d (nombre de décimales, défaut 1000)
 * -t (nombre de threads, défaut 1)
 * -s (afficher stats)
 * -q (pas afficher la valeur de pi, juste les stats)
 * -h (aide)
 * Version 1.0 - 06/01/2026 - Adrien Linuxtricks 
 * Version 1.1 - 06/01/2026 - Adrien Linuxtricks - Ajout stats
 * Version 2.0 - 10/01/2026 - Adrien Linuxtricks - Binary Splitting pour meilleurs perfs
 * Version 2.1 - 10/01/2026 - Adrien Linuxtricks - Ajout d'un mode quiet (ne pas afficher pi pour benchmarks par exemple)
 * Version 2.2 - 10/01/2026 - Adrien Linuxtricks - Ajout d'une aide
*/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <gmp.h>
#include <unistd.h>
#include <time.h>

// Structure pour stocker les résultats intermédiaires du binary splitting
// P, Q, T sont les composantes de la série de Chudnovsky décomposée
typedef struct {
    mpz_t P, Q, T;
} pqt_t;

// Structure pour passer les paramètres à chaque thread lors du calcul parallèle
typedef struct {
    int a, b;           // Intervalle [a, b) à calculer
    pqt_t result;       // Résultat P, Q, T pour cet intervalle
    int thread_id;      // Identifiant du thread
} bs_data_t;

// Affichage de l'aide
void print_help(const char *prog_name) {
    printf("Usage: %s [OPTIONS]\n\n", prog_name);
    printf("Calcul des décimales de PI via l'algorithme de Chudnovsky.\n\n");
    printf("OPTIONS:\n");
    printf("  -d NUM     Nombre de décimales à calculer (défaut: 1000)\n");
    printf("  -t NUM     Nombre de threads à utiliser (défaut: 1)\n");
    printf("  -s         Afficher les statistiques d'exécution (sur stderr)\n");
    printf("  -q         Mode silencieux (n'affiche pas la valeur de PI)\n");
    printf("  -h         Afficher cette aide et quitter\n\n");
    printf("EXEMPLES:\n");
    printf("  %s -d 5000                  Calculer 5000 décimales\n", prog_name);
    printf("  %s -d 10000 -t 4            Calculer 10000 décimales avec 4 threads\n", prog_name);
    printf("  %s -d 1000 -s               Afficher les statistiques\n", prog_name);
    printf("  %s -d 50000 -t 8 -s -q      Mode silencieux avec stats\n\n", prog_name);
}

// Fonction de binary splitting récursive pour calculer la série de Chudnovsky
// Cette méthode évite de recalculer les factorielles à chaque itération
// Divise l'intervalle [a, b) récursivement jusqu'aux cas de base
void binary_split(int a, int b, pqt_t *res) {
    if (b - a == 1) {
        // Cas de base : calcul d'un seul terme de la série
        mpz_t temp;
        mpz_init(temp);
        
        if (a == 0) {
            // Premier terme : initialisation P=1, Q=1
            mpz_set_ui(res->P, 1);
            mpz_set_ui(res->Q, 1);
        } else {
            // P(a,b) = -(6a-5)(2a-1)(6a-1)
            // Ce terme vient de la dérivation de (6k)!/(3k)!(k!)^3
            mpz_set_ui(res->P, 6 * a - 5);
            mpz_mul_ui(res->P, res->P, 2 * a - 1);
            mpz_mul_ui(res->P, res->P, 6 * a - 1);
            mpz_neg(res->P, res->P);
            
            // Q(a,b) = 10939058860032000 * a^3
            // Représente le dénominateur (640320)^(3k) factorisé
            mpz_set_ui(temp, a);
            mpz_pow_ui(temp, temp, 3);
            mpz_mul_ui(res->Q, temp, 640320);
            mpz_mul_ui(res->Q, res->Q, 640320);
            mpz_mul_ui(res->Q, res->Q, 640320);
            mpz_divexact_ui(res->Q, res->Q, 24);
        }
        
        // T(a,b) = P(a,b) * (13591409 + 545140134*a)
        // Le numérateur de la série avec les constantes de Chudnovsky
        mpz_set_ui(temp, 545140134);
        mpz_mul_ui(temp, temp, a);
        mpz_add_ui(temp, temp, 13591409);
        mpz_mul(res->T, res->P, temp);
        
        mpz_clear(temp);
    } else {
        // Cas récursif : diviser l'intervalle en deux parties
        int m = (a + b) / 2;
        
        pqt_t left, right;
        mpz_inits(left.P, left.Q, left.T, NULL);
        mpz_inits(right.P, right.Q, right.T, NULL);
        
        // Calcul récursif des deux moitiés
        binary_split(a, m, &left);
        binary_split(m, b, &right);
        
        // Combinaison des résultats gauche et droite selon les formules de binary splitting
        mpz_t temp;
        mpz_init(temp);
        
        // P(a,b) = P(a,m) * P(m,b)
        mpz_mul(res->P, left.P, right.P);
        
        // Q(a,b) = Q(a,m) * Q(m,b)
        mpz_mul(res->Q, left.Q, right.Q);
        
        // T(a,b) = Q(m,b) * T(a,m) + P(a,m) * T(m,b)
        // Cette formule permet de combiner les sommes partielles
        mpz_mul(res->T, right.Q, left.T);
        mpz_mul(temp, left.P, right.T);
        mpz_add(res->T, res->T, temp);
        
        // Libération de la mémoire des calculs intermédiaires
        mpz_clear(temp);
        mpz_clears(left.P, left.Q, left.T, NULL);
        mpz_clears(right.P, right.Q, right.T, NULL);
    }
}

// Fonction exécutée par chaque thread pour calculer une partie de la série via binary splitting
void *compute_bs_thread(void *arg) {
    bs_data_t *data = (bs_data_t *)arg;
    binary_split(data->a, data->b, &data->result);
    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    int decimals = 1000; // Par défaut ou repris si option -d (voir dessous)
    int num_threads = 1; // Par défaut ou repris si option -t (voir dessous)
    int show_stats = 0;  // Par défaut pas de stats sauf si actif via option -s (voir dessous)
    int quiet_mode = 0;  // Par défaut afficher PI sauf si actif via option -q (voir dessous)
    int opt;

    // Analyse des options : -d pour décimales, -t pour threads, -s pour stats, -q pour quiet, -h pour aide
    while ((opt = getopt(argc, argv, "d:t:sqh")) != -1) {
        switch (opt) {
            case 'd':
                decimals = atoi(optarg);
                break;
            case 't':
                num_threads = atoi(optarg);
                if (num_threads < 1) num_threads = 1;
                break;
            case 's':
                show_stats = 1;
                break;
            case 'q':
                quiet_mode = 1;
                break;
            case 'h':
                print_help(argv[0]);
                return 0;
            default:
                print_help(argv[0]);
                return 1;
        }
    }

    // Début du chronométrage
    struct timespec start_time, end_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    // Précision GMP augmentée pour la stabilité (environ 4 bits par décimale)
    mpf_set_default_prec((decimals + 100) * 4);

    // Nombre d'itérations nécessaires pour Chudnovsky (environ 14 décimales par itération)
    int iterations = decimals / 14 + 10;

    if (num_threads == 1) {
        // Version mono-thread : calcul direct via binary splitting
        pqt_t result;
        mpz_inits(result.P, result.Q, result.T, NULL);
        
        binary_split(0, iterations, &result);
        
        // Conversion en nombres flottants pour le calcul final de PI
        mpf_t pi, fP, fQ, fT, sqrt_C;
        mpf_inits(pi, fP, fQ, fT, sqrt_C, NULL);
        
        mpf_set_z(fP, result.P);
        mpf_set_z(fQ, result.Q);
        mpf_set_z(fT, result.T);
        
        // Calcul final de PI basé sur la formule : PI = (426880 * sqrt(10005)) * Q / T
        mpf_sqrt_ui(sqrt_C, 10005);
        mpf_mul_ui(sqrt_C, sqrt_C, 426880);
        mpf_mul(pi, sqrt_C, fQ);
        mpf_div(pi, pi, fT);
        
        if (!quiet_mode) {
            gmp_printf("%.*Ff\n", decimals, pi);
        }
        
        // Libération de la mémoire
        mpz_clears(result.P, result.Q, result.T, NULL);
        mpf_clears(pi, fP, fQ, fT, sqrt_C, NULL);
    } else {
        // Version multi-thread : répartition du travail entre threads
        int chunk = iterations / num_threads;
        
        // Allocation dynamique selon le nombre de threads
        pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
        bs_data_t *data = malloc(num_threads * sizeof(bs_data_t));
        
        // Création des threads et répartition du travail
        for (int i = 0; i < num_threads; i++) {
            data[i].a = i * chunk;
            data[i].b = (i == num_threads - 1) ? iterations : (i + 1) * chunk;
            data[i].thread_id = i;
            mpz_inits(data[i].result.P, data[i].result.Q, data[i].result.T, NULL);
            pthread_create(&threads[i], NULL, compute_bs_thread, &data[i]);
        }
        
        // Attente des threads et combinaison des résultats
        pqt_t final;
        mpz_inits(final.P, final.Q, final.T, NULL);
        
        // Initialisation avec le premier résultat
        pthread_join(threads[0], NULL);
        mpz_set(final.P, data[0].result.P);
        mpz_set(final.Q, data[0].result.Q);
        mpz_set(final.T, data[0].result.T);
        
        // Agrégation des résultats des autres threads
        for (int i = 1; i < num_threads; i++) {
            pthread_join(threads[i], NULL);
            
            mpz_t temp;
            mpz_init(temp);
            
            // Combiner les résultats selon les formules de binary splitting
            mpz_mul(temp, data[i].result.Q, final.T);
            mpz_mul(final.T, final.P, data[i].result.T);
            mpz_add(final.T, final.T, temp);
            
            mpz_mul(final.P, final.P, data[i].result.P);
            mpz_mul(final.Q, final.Q, data[i].result.Q);
            
            mpz_clear(temp);
            mpz_clears(data[i].result.P, data[i].result.Q, data[i].result.T, NULL);
        }
        
        // Calcul final de PI basé sur la somme obtenue
        mpf_t pi, fQ, fT, sqrt_C;
        mpf_inits(pi, fQ, fT, sqrt_C, NULL);
        
        mpf_set_z(fQ, final.Q);
        mpf_set_z(fT, final.T);
        
        mpf_sqrt_ui(sqrt_C, 10005);
        mpf_mul_ui(sqrt_C, sqrt_C, 426880);
        mpf_mul(pi, sqrt_C, fQ);
        mpf_div(pi, pi, fT);
        
        if (!quiet_mode) {
            gmp_printf("%.*Ff\n", decimals, pi);
        }
        
        // Nettoyage de la mémoire
        mpz_clears(final.P, final.Q, final.T, NULL);
        mpf_clears(pi, fQ, fT, sqrt_C, NULL);
        mpz_clears(data[0].result.P, data[0].result.Q, data[0].result.T, NULL);
        free(threads);
        free(data);
    }

    // Fin du chronométrage
    clock_gettime(CLOCK_MONOTONIC, &end_time);

    // Affichage des statistiques si l'option -s est activée mais pas sur stdout si redirection dans fichier ou /dev/null
    if (show_stats) {
        double elapsed = (end_time.tv_sec - start_time.tv_sec) + 
                        (end_time.tv_nsec - start_time.tv_nsec) / 1e9;
        double decimals_per_second = decimals / elapsed;

        fprintf(stderr, "======= Stats =======\n");
        fprintf(stderr, "Time      : %.3f s\n", elapsed);
        fprintf(stderr, "Threads   : %d\n", num_threads);
        fprintf(stderr, "Decimals  : %d\n", decimals);
        fprintf(stderr, "Dec / sec : %.0f\n", decimals_per_second);
    }

    return 0;
}
