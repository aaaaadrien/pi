# pi
Programme C pour calculer PI

## Dependances
```
dnf install gmp-devel
```

## Compilation
```
gcc -O3 -o pi pi.c -lgmp -lpthread -lm
```

## Usage
```
Usage: pi [OPTIONS]

Calcul des décimales de PI via l'algorithme de Chudnovsky.

OPTIONS:
  -d NUM     Nombre de décimales à calculer (défaut: 1000)
  -t NUM     Nombre de threads à utiliser (défaut: 1)
  -s         Afficher les statistiques d'exécution (sur stderr)
  -q         Mode silencieux (n'affiche pas la valeur de PI)
  -h         Afficher cette aide et quitter

EXEMPLES:
  ./pi -d 5000                  Calculer 5000 décimales
  ./pi -d 10000 -t 4            Calculer 10000 décimales avec 4 threads
  ./pi -d 1000 -s               Afficher les statistiques
  ./pi -d 50000 -t 8 -s -q      Mode silencieux avec stats
```
