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

Pour exécuter notre code, ouvrez le fichier run.sh, sélectionez l'instance sur laquelle exécuter le code (tout les fichiers sont dans le dossier data/processed), sélectionnez une des méthode disponible, choisissez la time limit (en secondes) pour les solveurs ainsi que le niveau de verbose (entre 0 et 2). Enfin, lancez l'exécution du fichier run.sh depuis l'invite de commande.

## Structure du projet

Le projet est réalisé en C++, seule une petite partie de post-traitement est réalisée en python.

Pour parser les instances, nous avons créé une classe _Instance_, décrite dans le parser.h. Cet objet contient toutes les informations requises concernant l'instance pour l'éxécution des différents algorithmes par la suite. On peut également stocker dans l'objet une solution. La classe dispose également de méthodes permettant de calculer les scores et contraintes statiques et robustes correspondant à la solution stockée dans l'objet.

Pour la résolution, nous avons créé une classe mère _SolveMethode_ dont héritent les classes gérant la résolution par les différentes méthodes. La classe sert de support pour contenir les méthodes utilisées par plusieurs méthodes de résolution. 
La structure du code pour les méthodes statique, dualisé, plans coupants et branch and cut sont très similaires. Pour l'heuristique, en plus de la classe permettant la résolution de l'instance par notre algorithme, nous avons ajouté deux classes. L'une permet de stocker efficacement les information nécessaires lors de l'exploration du graphe par A*, et l'autre garde en mémoire les solutions obtenue lors de la recherche de la pondération optimale.

Les résultats sont stockés dans le dossier results sous la forme de CSV.