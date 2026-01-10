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
pi [-d decimals] [-t threads]
```
