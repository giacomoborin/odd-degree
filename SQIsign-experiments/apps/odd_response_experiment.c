// SPDX-License-Identifier: Apache-2.0

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bench_test_arguments.h>
#include <encoded_sizes.h>
#include <rng.h>
#include <signature.h>
#if defined(TARGET_BIG_ENDIAN)
#include <tutil.h>
#endif

#define STR_(x) #x
#define STR(x) STR_(x)

static int
parse_u64_option(const char *arg, const char *name, uint64_t *out)
{
    size_t len = strlen(name);
    if (strncmp(arg, name, len) != 0 || arg[len] != '=')
        return 0;

    errno = 0;
    char *end = NULL;
    unsigned long long value = strtoull(arg + len + 1, &end, 10);
    if (errno != 0 || end == arg + len + 1 || *end != '\0') {
        fprintf(stderr, "Invalid integer option: %s\n", arg);
        exit(1);
    }

    *out = (uint64_t)value;
    return 1;
}

static void
print_header(void)
{
    printf("pair,total_responses,odd_responses,odd_responses/total_responses,complete,candidates");
    printf("\n");
}

static void
print_stats(uint64_t index, const response_isogeny_stats_t *stats)
{
    double odd_ratio = stats->responses == 0 ? 0.0 : ((double)stats->odd_degree) / (double)stats->responses;

    printf("%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%.10f,%u,%" PRIu64,
           index,
           stats->responses,
           stats->odd_degree,
           odd_ratio,
           stats->complete,
           stats->candidates);
    printf("\n");
}

static void
usage(const char *argv0)
{
    printf("Usage: %s [--iterations=<n>] [--response-length=<e>] [--max-responses=<n>] [--max-candidates=<n>] [--seed=<seed>]\n",
           argv0);
    printf("  --iterations       commitment/challenge pairs to sample under one key (default: 1)\n");
    printf("  --max-candidates   parallelogram candidates/attempts per pair; 0 means unlimited (default: 1000000)\n");
    printf("  --max-responses    random valid responses to sample per pair; 0 means exact enumeration (default: 0)\n");
    printf("  --response-length  response exponent used in the enumeration bound (default: %d)\n", SQIsign_response_length);
    printf("  --seed             deterministic 12-word seed, same format as the test binaries\n");
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    uint64_t iterations = 1;
    uint64_t max_candidates = 1000000;
    uint64_t max_responses = 0;
    uint64_t response_length_arg = SQIsign_response_length;
    int seed_set = 0;
    int help = 0;

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }

        if (parse_u64_option(argv[i], "--iterations", &iterations))
            continue;

        if (parse_u64_option(argv[i], "--max-candidates", &max_candidates))
            continue;

        if (parse_u64_option(argv[i], "--max-responses", &max_responses))
            continue;

        if (parse_u64_option(argv[i], "--response-length", &response_length_arg))
            continue;

        fprintf(stderr, "Unknown option: %s\n", argv[i]);
        return 1;
    }

    if (help || iterations == 0) {
        usage(argv[0]);
        return help ? 0 : 1;
    }
    if (response_length_arg > UINT32_MAX) {
        fprintf(stderr, "--response-length is too large\n");
        return 1;
    }
    uint32_t response_length = (uint32_t)response_length_arg;

    if (!seed_set) {
        randombytes_select((unsigned char *)seed, sizeof(seed));
    }

    print_seed(seed);

#if defined(TARGET_BIG_ENDIAN)
    for (int i = 0; i < 12; i++) {
        seed[i] = BSWAP32(seed[i]);
    }
#endif

    randombytes_init((unsigned char *)seed, NULL, 256);

    public_key_t pk;
    secret_key_t sk;
    public_key_init(&pk);
    secret_key_init(&sk);

    if (!protocols_keygen(&pk, &sk)) {
        fprintf(stderr, "Key generation failed\n");
        secret_key_finalize(&sk);
        public_key_finalize(&pk);
        return 1;
    }

    const unsigned char message[32] = { 0 };
    response_isogeny_stats_t aggregate = { 0 };
    aggregate.complete = 1;

    printf("# variant=%s response_length=%" PRIu32 " max_candidates=%" PRIu64
           " max_responses=%" PRIu64 "\n",
           STR(SQISIGN_VARIANT),
           response_length,
           max_candidates,
           max_responses);
    print_header();
    fflush(stdout);

    for (uint64_t i = 0; i < iterations; i++) {
        response_isogeny_stats_t stats;
        if (!protocols_response_isogeny_stats(
                &stats,
                &pk,
                &sk,
                message,
                sizeof(message),
                response_length,
                max_candidates,
                max_responses)) {
            fprintf(stderr, "Experiment failed for pair %" PRIu64 "\n", i);
            secret_key_finalize(&sk);
            public_key_finalize(&pk);
            return 1;
        }

        print_stats(i, &stats);
        fflush(stdout);

        aggregate.candidates += stats.candidates;
        aggregate.responses += stats.responses;
        aggregate.odd_degree += stats.odd_degree;
        aggregate.even_degree += stats.even_degree;
        aggregate.invalid += stats.invalid;
        aggregate.complete &= stats.complete;
    }

    double odd_ratio = aggregate.responses == 0 ? 0.0 : ((double)aggregate.odd_degree) / (double)aggregate.responses;
    printf("# aggregate candidates=%" PRIu64 " total_responses=%" PRIu64 " odd_responses=%" PRIu64
           " even_responses=%" PRIu64 " invalid=%" PRIu64 " odd_ratio=%.10f complete=%u\n",
           aggregate.candidates,
           aggregate.responses,
           aggregate.odd_degree,
           aggregate.even_degree,
           aggregate.invalid,
           odd_ratio,
           aggregate.complete);
    fflush(stdout);

    secret_key_finalize(&sk);
    public_key_finalize(&pk);

    return 0;
}
