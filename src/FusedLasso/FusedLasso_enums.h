#pragma once

/**
 * @file FusedLasso_enums.h
 * @brief Énumérés centralisés pour les types de régression et de pénalité
 */

/**
 * @enum regEnum
 * @brief Types de régression supportés
 */
enum class regEnum {
    GAUSSIAN,  ///< Régression linéaire (Gaussian)
    BINOMIAL   ///< Régression logistique (Binomial)
};

/**
 * @enum class penEnum
 * @brief Types de pénalité supportés
 */
enum class penEnum {
    L1,        ///< Pénalité L1 (norme L1)
    Huber,     ///< Pénalité de Huber (pénalité robuste)
    L2         ///< Pénalité L2 (norme L2)
};

/**
 * @enum FusionStrategy
 * @brief Niveau maximal de fusion autorisé durant l'algorithme
 */
enum class FusionStrategy {
    None              = -1,  ///< Aucune fusion (beta figé dans les groupes courants)
    EqualOnly         =  0,  ///< Fusionne uniquement les groupes égaux
    EqualSplitInactive=  1,  ///< Ajoute le split des inactifs
    Full              =  2   ///< Stratégie complète : égaux + split inactifs + split actifs
};

