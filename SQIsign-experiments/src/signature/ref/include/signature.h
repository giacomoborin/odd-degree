/** @file
 *
 * @brief The key generation and signature protocols
 */

#ifndef SIGNATURE_H
#define SIGNATURE_H

#include <stddef.h>
#include <stdint.h>

#include <sqisign_namespace.h>
#include <ec.h>
#include <quaternion.h>
#include <verification.h>

/** @defgroup signature SQIsignHD key generation and signature protocols
 * @{
 */
/** @defgroup signature_t Types for SQIsignHD key generation and signature protocols
 * @{
 */

/** @brief Type for the secret keys
 *
 * @typedef secret_key_t
 *
 * @struct secret_key
 *
 */
typedef struct secret_key
{
    ec_curve_t curve; /// the public curve, but with little precomputations
    quat_left_ideal_t secret_ideal;
    ibz_mat_2x2_t mat_BAcan_to_BA0_two; // mat_BA0_to_BAcan*BA0 = BAcan, where BAcan is the
                                        // canonical basis of EA[2^e], and BA0 the image of the
                                        // basis of E0[2^e] through the secret isogeny
    ec_basis_t canonical_basis;         // the canonical basis of the public key curve
} secret_key_t;

/** @brief Statistics for enumerating response isogenies for one commitment/challenge pair */
typedef struct response_isogeny_stats
{
    uint64_t candidate_space;
    uint64_t candidates;
    uint64_t responses;
    uint64_t odd_degree;
    uint64_t even_degree;
    uint64_t invalid;
    uint8_t candidate_space_overflow;
    uint8_t complete;
} response_isogeny_stats_t;

/** @}
 */

/*************************** Functions *****************************/

void secret_key_init(secret_key_t *sk);
void secret_key_finalize(secret_key_t *sk);

/**
 * @brief Key generation
 *
 * @param pk Output: will contain the public key
 * @param sk Output: will contain the secret key
 * @returns 1 if success, 0 otherwise
 */
int protocols_keygen(public_key_t *pk, secret_key_t *sk);

/**
 * @brief Signature computation
 *
 * @param sig Output: will contain the signature
 * @param sk secret key
 * @param pk public key
 * @param m message
 * @param l size
 * @returns 1 if success, 0 otherwise
 */
int protocols_sign(signature_t *sig, const public_key_t *pk, secret_key_t *sk, const unsigned char *m, size_t l);

/**
 * @brief Enumerate response isogeny candidates for one freshly sampled commitment/challenge pair.
 *
 * The function samples a commitment using the regular SQIsign commitment code,
 * hashes the supplied message to the corresponding challenge, enumerates the
 * same response lattice ball used by signing, with the supplied response
 * exponent, and records how many response isogenies have odd degree. Set
 * max_candidates and/or max_responses to 0 for an unbounded exact enumeration.
 *
 * @param stats Output statistics
 * @param pk public key
 * @param sk secret key
 * @param m message
 * @param l size
 * @param response_length exponent e_rsp used in the response norm bound
 * @param max_candidates maximum parallelogram candidates to inspect, or 0 for no cap
 * @param max_responses maximum valid responses to inspect, or 0 for no cap
 * @returns 1 if the commitment/challenge pair was generated and enumeration ran, 0 otherwise
 */
int protocols_response_isogeny_stats(response_isogeny_stats_t *stats,
                                     const public_key_t *pk,
                                     secret_key_t *sk,
                                     const unsigned char *m,
                                     size_t l,
                                     uint32_t response_length,
                                     uint64_t max_candidates,
                                     uint64_t max_responses);

/*************************** Encoding *****************************/

/** @defgroup encoding Encoding and decoding functions
 * @{
 */

/**
 * @brief Encodes a secret key as a byte array
 *
 * @param enc : Byte array to encode the secret key (including public key) in
 * @param sk : Secret key to encode
 * @param pk : Public key to encode
 */
void secret_key_to_bytes(unsigned char *enc, const secret_key_t *sk, const public_key_t *pk);

/**
 * @brief Decodes a secret key (and public key) from a byte array
 *
 * @param sk : Structure to decode the secret key in
 * @param pk : Structure to decode the public key in
 * @param enc : Byte array to decode
 */
void secret_key_from_bytes(secret_key_t *sk, public_key_t *pk, const unsigned char *enc);

/** @}
 */

/** @}
 */

#endif
