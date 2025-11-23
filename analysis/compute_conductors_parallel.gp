/*
 * Compute log(conductors) in parallel using PARI/GP
 *
 * Reads data/ainvs_of_big_starter_dataset.txt and computes
 * the log of the conductor for each elliptic curve using parallel processing.
 *
 * Input format: rational,a1,a2,a3,a4,a6
 * Output format: rational,log_conductor
 *
 * Uses a timeout mechanism to handle curves with very large conductors
 * that are slow to factor. If computation exceeds CONDUCTOR_TIMEOUT seconds,
 * the log conductor is set to BIG_LOG_CONDUCTOR (default: 10000).
 *
 * Configuration parameters (edit lines 18-20):
 *   CONDUCTOR_TIMEOUT: timeout in seconds (default: 3)
 *   BIG_LOG_CONDUCTOR: default value for timed-out curves (default: 10000)
 *   MAX_CURVES: limit processing to first N curves; set to 0 to process all (default: 10)
 */

default(threadsizemax, 100000000)

/* Global configuration parameters */
CONDUCTOR_TIMEOUT = 3;        /* Timeout in seconds - adjust as needed */
BIG_LOG_CONDUCTOR = 10000;    /* Default value for timeout cases */
MAX_CURVES = 0;              /* Max curves to process (0 = process all) - set to small number for testing */

/* Function to compute log(conductor) from ainvs - must be self-contained for parapply */
compute_conductor(ainvs_data) =
{
    my(rational, ainvs, E, red, conductor, log_conductor, timeout_seconds, big_value);

    /* Local copies of global constants (needed for parapply) */
    timeout_seconds = 3;        /* CONDUCTOR_TIMEOUT */
    big_value = 10000;          /* BIG_LOG_CONDUCTOR */

    /* ainvs_data is a vector [rational, a1, a2, a3, a4, a6] */
    rational = ainvs_data[1];
    ainvs = [ainvs_data[2], ainvs_data[3], ainvs_data[4], ainvs_data[5], ainvs_data[6]];

    /* Try to compute conductor with timeout */
    iferr(alarm(timeout_seconds); E = ellinit(ainvs); red = ellglobalred(E); conductor = red[1]; log_conductor = log(conductor); alarm(0), err, alarm(0); print("Timeout on ", rational, ", using default value ", big_value); log_conductor = big_value);

    /* Return [rational, log_conductor] as a vector */
    [rational, log_conductor]
}

print("Reading data/ainvs_of_big_starter_dataset.txt...");

/* Read the entire file */
lines = readstr("../data/ainvs_of_big_starter_dataset.txt");
total_curves = length(lines);

/* Limit number of curves if MAX_CURVES > 0 */
if(MAX_CURVES > 0 && MAX_CURVES < total_curves, n = MAX_CURVES; print("Found ", total_curves, " curves, processing first ", n, " for testing"), n = total_curves; print("Found ", n, " curves, processing all"));

print("Parsing input data...");

/* Parse each line into [rational, a1, a2, a3, a4, a6] */
parsed_data = vector(n);
for(j=1, n, my(parts = strsplit(lines[j], ",")); parsed_data[j] = [parts[1], eval(parts[2]), eval(parts[3]), eval(parts[4]), eval(parts[5]), eval(parts[6])]);

print("Parsed ", n, " lines");
print("Computing log(conductors) in parallel (using ", default(nbthreads), " threads)...");
print("");

/* Apply compute_conductor in parallel to all curves */
results = parapply(compute_conductor, parsed_data);

print("");
print("Writing results to data/log_conductors_of_big_starter_dataset.txt...");

/* Write results to file */
write1("../data/log_conductors_of_big_starter_dataset.txt", "");  /* Clear file */

for(j=1, n, write("../data/log_conductors_of_big_starter_dataset.txt", results[j][1], ",", results[j][2]));

print("");
print("Done! Processed ", n, " curves.");

/* Optional: show some statistics */
log_conductors = vector(n, k, results[k][2]);
print("");
print("Statistics:");
print("  Minimum log(conductor): ", vecmin(log_conductors));
print("  Maximum log(conductor): ", vecmax(log_conductors));

quit;
