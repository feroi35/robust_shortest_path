# Robust quickest path

Projet du MPRO
2023-2024
De Zacharie Alès et Daniel Porumbel

Benoit Duval et Raphaël Taisant

Le but de ce projet est l'étude d'un problème de plus court chemin robuste. On résoudra ce problème en utilisant:

- Un algorithme de plans coupants
- Un algorithme de Branch-and-Cut (via un LazyCallback)
- Une méthode de dualisation
- Une heuristique

Pour tester les résolutions, on utilisera les instances de la 9e édition de DIMACS avec des villes américaines avec entre 20 et 3500 sommets.

http://www.diag.uniroma1.it/~challenge9/
