#pragma once

#include "FusedLasso_utils.h"
#include <vector>

/**
 * @file FusedLasso_vectorops.h
 * @brief Utilitaires pour opérations vectorielles communes
 * 
 * Centralise les opérations vectorielles fréquemment utilisées
 * pour éviter la duplication de code entre les différentes classes.
 */

namespace VectorOps {

    /**
     * @brief Multiplie un vecteur par un scalaire en place
     * @param vec Vecteur à modifier
     * @param scalar Valeur scalaire à multiplier
     */
    inline void scalarMultiply(vector<double>& vec, double scalar) {
        for (auto& v : vec) v *= scalar ;
    }

    /**
     * @brief Multiplie les éléments de vec par les poids correspondants en place
     * @param vec Vecteur à modifier
     * @param weights Vecteur de poids (doit avoir la même taille que vec)
     */
    inline void weightedMultiply(vector<double>& vec, const vector<double>& weights) {
        if (vec.size() != weights.size()) {
            Rcpp::stop("VectorOps::weightedMultiply: Tailles incompatibles\n");
        }
        for (size_t i = 0; i < vec.size(); ++i) vec[i] *= weights[i] ;
    }

    /**
     * @brief Ajoute un vecteur à un autre en place: result = result + scalar * other
     * @param result Vecteur résultat (modifié en place)
     * @param other Vecteur à ajouter
     * @param scalar Facteur multiplicatif (par défaut 1.0)
     */
    inline void addScaledVector(vector<double>& result, const vector<double>& other, double scalar = 1.0) {
        if (result.size() != other.size()) {
            Rcpp::stop("VectorOps::addScaledVector: Tailles incompatibles\n");
        }
        for (size_t i = 0; i < result.size(); ++i) result[i] += scalar * other[i] ;
    }

    /**
     * @brief Calcule le produit scalaire de deux vecteurs
     * @param vec1 Premier vecteur
     * @param vec2 Deuxième vecteur
     * @return Produit scalaire
     */
    inline double dotProduct(const vector<double>& vec1, const vector<double>& vec2) {
        if (vec1.size() != vec2.size()) {
            Rcpp::stop("VectorOps::dotProduct: Tailles incompatibles\n");
        }
        double result = 0.0;
        for (size_t i = 0; i < vec1.size(); ++i) result += vec1[i] * vec2[i] ;
        return result;
    }

    /**
     * @brief Calcule le produit scalaire pondéré: sum(vec1[i] * vec2[i] * weights[i])
     * @param vec1 Premier vecteur
     * @param vec2 Deuxième vecteur
     * @param weights Vecteur de poids
     * @return Produit scalaire pondéré
     */
    inline double weightedDotProduct(const vector<double>& vec1, const vector<double>& vec2, const vector<double>& weights) {
        if (vec1.size() != vec2.size() || vec1.size() != weights.size()) {
            Rcpp::stop("VectorOps::weightedDotProduct: Tailles incompatibles\n");
        }
        double result = 0.0;
        for (size_t i = 0; i < vec1.size(); ++i) result += vec1[i] * vec2[i] * weights[i] ;
        return result;
    }

    /**
     * @brief Initialise un vecteur avec une valeur constante
     * @param vec Vecteur à initialiser
     * @param value Valeur à affecter à chaque élément
     */
    inline void fill(vector<double>& vec, double value) {
        for (auto& v : vec) v = value ;
    }

    /**
     * @brief Définit les éléments d'un vecteur de manière élément par élément
     * avec la fonction f: vec[i] = f(vec[i])
     * @param vec Vecteur à transformer
     * @param f Fonction de transformation
     */
    template<typename F>
    inline void apply(vector<double>& vec, F f) {
        for (auto& v : vec) v = f(v) ;
    }

    /**
     * @brief Compte le nombre de zéros dans un vecteur
     * @param vec Vecteur à analyser
     * @return Nombre d'éléments égaux à zéro
     */
    inline int countZeros(const vector<double>& vec) {
        int count = 0;
        for (const auto& v : vec) if (v == 0.0) count++ ;
        return count;
    }

    /**
     * @brief Compte les éléments non zéro dans un vecteur
     * @param vec Vecteur à analyser
     * @return Nombre d'éléments différents de zéro
     */
    inline int countNonZeros(const vector<double>& vec) {
        return vec.size() - countZeros(vec);
    }

    /**
     * @brief Écrête les valeurs du vecteur entre min et max
     * @param vec Vecteur à modifier
     * @param minVal Valeur minimum
     * @param maxVal Valeur maximum
     */
    inline void clamp(vector<double>& vec, double minVal, double maxVal) {
        for (auto& v : vec) { if (v < minVal) v = minVal; else if (v > maxVal) v = maxVal; }
    }

}

