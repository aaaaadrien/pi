/* * Programme de calcul des décimales de PI via l'algorithme de Chudnovsky.
 * Utilise la bibliothèque GMP (assez précis) et Pthreads pour le multithread.
 * Options : -d (nombre de décimales, défaut 1000) et -t (nombre de threads, défaut 1).
 * Version 1.0 - 06/01/2026 - Adrien Linuxtricks 
*/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <gmp.h>
#include <unistd.h>

// Structure pour passer les paramètres à chaque thread
typedef struct {
    int start;
    int end;
    mpf_t sum;
} thread_data_t;

// Fonction exécutée par chaque thread pour calculer une partie de la série
void *compute_partial_sum(void *arg) {
    thread_data_t *data = (thread_data_t *)arg;

    mpf_init(data->sum);
    mpf_set_ui(data->sum, 0);

    mpz_t num, den, fact6k, fact3k, factk, pow;
    mpf_t term, fnum, fden;

    mpz_inits(num, den, fact6k, fact3k, factk, pow, NULL);
    mpf_inits(term, fnum, fden, NULL);

    for (int k = data->start; k < data->end; k++) {
        // Formule de Chudnovsky : calcul des factorielles https://fr.wikipedia.org/wiki/Algorithme_de_Chudnovski
        mpz_fac_ui(fact6k, 6 * k);
        mpz_fac_ui(fact3k, 3 * k);
        mpz_fac_ui(factk, k);

        // num = (6k)! * (13591409 + 545140134k)
        mpz_mul_ui(num, fact6k, 13591409 + 545140134UL * k);

        // den = (3k)! * (k!)^3 * (640320)^(3k)
        mpz_pow_ui(pow, factk, 3);
        mpz_mul(den, fact3k, pow);
        mpz_ui_pow_ui(pow, 640320, 3 * k);
        mpz_mul(den, den, pow);

        mpf_set_z(fnum, num);
        mpf_set_z(fden, den);
        mpf_div(term, fnum, fden);

        // Alternance du signe (-1)^k
        if (k % 2) mpf_neg(term, term);

        mpf_add(data->sum, data->sum, term);
    }

    // Libération de la mémoire locale au thread
    mpz_clears(num, den, fact6k, fact3k, factk, pow, NULL);
    mpf_clears(term, fnum, fden, NULL);

    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    int decimals = 1000; // Par défaut ou repris si opton -d (voir dessous)
    int num_threads = 1; // Par défaut ou repris si opton -t (voir dessous)
    int opt;

    // Analyse des options : -d pour décimales, -t pour threads
    while ((opt = getopt(argc, argv, "d:t:")) != -1) {
        switch (opt) {
            case 'd':
                decimals = atoi(optarg);
                break;
            case 't':
                num_threads = atoi(optarg);
                if (num_threads < 1) num_threads = 1;
                break;
        }
    }

    // Précision GMP (environ 4 bits par décimale)
    mpf_set_default_prec(decimals * 4);

    int iterations = decimals / 14 + 2;
    int chunk = iterations / num_threads;

    // Allocation dynamique selon le nombre de threads
    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    thread_data_t *data = malloc(num_threads * sizeof(thread_data_t));

    // Création des threads et répartition du travail
    for (int i = 0; i < num_threads; i++) {
        data[i].start = i * chunk;
        data[i].end = (i == num_threads - 1) ? iterations : (i + 1) * chunk;
        pthread_create(&threads[i], NULL, compute_partial_sum, &data[i]);
    }

    mpf_t sum;
    mpf_init(sum);
    mpf_set_ui(sum, 0);

    // Attente des threads et agrégation des résultats
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
        mpf_add(sum, sum, data[i].sum);
        mpf_clear(data[i].sum);
    }

    // Calcul final de PI basé sur la somme obtenue
    mpf_t pi, C, sqrt10005;
    mpf_inits(pi, C, sqrt10005, NULL);

    mpf_sqrt_ui(sqrt10005, 10005);
    mpf_mul_ui(C, sqrt10005, 426880);
    mpf_div(pi, C, sum);

    gmp_printf("%.*Ff\n", decimals, pi);

    // Nettoyage de la mémoire
    mpf_clears(sum, pi, C, sqrt10005, NULL);
    free(threads);
    free(data);

    return 0;
}
